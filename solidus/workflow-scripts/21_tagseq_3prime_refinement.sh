#!/usr/bin/env bash
# Reconstruct the TagSeq-supported 3-prime annotation refinement used in 2026.
#
# Historical software versions:
#   cutadapt 5.2; STAR 2.7.2a; samtools 1.21; bedtools 2.31.1.
#
# The TagSeq libraries are forward/same-strand (98.34% in the C9 pilot).
# Only reads containing a detected poly(A) tract are used. The workflow changes
# transcript and terminal-exon boundaries only; CDS coordinates are unchanged.
set -euo pipefail

TAG=${TAG:-/work/shared/ssol-annotation/Bolnick_worm_Tagseq}
ROOT=${ROOT:-/work/cyu/annotation/ssol-annotation}
OUT=${OUT:-${ROOT}/results/tagseq_3prime}
POLY=${POLY:-${OUT}/polyA_sites}
GENOME=${GENOME:-/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa}
INPUT_GTF=${INPUT_GTF:-${ROOT}/results/tsa_alignment/ssol_braker_tsa_utr_plus_187.gtf}
FINAL_GTF=${FINAL_GTF:-${ROOT}/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.gtf}
CHECKED_GTF=${CHECKED_GTF:-${ROOT}/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.checked.gtf}
STAR_BIN=${STAR_BIN:-/work/cyu/scRNA/cellranger-10.0.0/lib/bin/STAR}
STARREF=${STARREF:-/work/cyu/scRNA/reference_test/ssol_braker_tsa_utr_plus_187/star}
THREADS=${THREADS:-24}

for cmd in cutadapt samtools bedtools awk sort gzip zcat md5sum; do
    command -v "${cmd}" >/dev/null || {
        echo "ERROR: ${cmd} is not in PATH" >&2
        exit 1
    }
done
# ARGIND in the annotation update requires GNU awk, not BSD awk or mawk.
awk --version 2>/dev/null | grep -q 'GNU Awk' || {
    echo 'ERROR: GNU awk must be available as awk; no jobs launched.' >&2
    exit 1
}
for file in "${FINAL_GTF}" "${CHECKED_GTF}"; do
    [[ ! -e "${file}" && ! -L "${file}" ]] || {
        echo "ERROR: output exists; use a new output path: ${file}" >&2
        exit 1
    }
done
[[ -x "${STAR_BIN}" ]] || { echo "ERROR: missing STAR: ${STAR_BIN}" >&2; exit 1; }
for file in "${GENOME}" "${INPUT_GTF}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done
[[ -d "${TAG}" ]] || { echo "ERROR: missing TagSeq directory: ${TAG}" >&2; exit 1; }
[[ -d "${STARREF}" ]] || { echo "ERROR: missing STAR index: ${STARREF}" >&2; exit 1; }

mkdir -p "${OUT}/polyA_reads" "${OUT}/logs" "${OUT}/star_align" "${POLY}"
ulimit -n 65535 2>/dev/null || true

# Samples have two lanes. Lane FASTQs are processed separately and then merged.
find "${TAG}" -maxdepth 1 -type f -name '*_L001_R1_001.fastq.gz' -print |
sed -E 's#^.*/##; s/_S[0-9]+_L001_R1_001.fastq.gz$//' |
sort -u > "${OUT}/sample_ids.txt"

echo "Samples: $(wc -l < "${OUT}/sample_ids.txt")"
[[ -s "${OUT}/sample_ids.txt" ]] || { echo 'ERROR: no matching TagSeq samples' >&2; exit 1; }

while read -r sample; do
    lane_outputs=()
    for raw in "${TAG}"/"${sample}"_S*_L00?_R1_001.fastq.gz; do
        [[ -e "${raw}" ]] || continue
        base=$(basename "${raw}" .fastq.gz)
        lane_out="${OUT}/polyA_reads/${base}.polyA.fastq.gz"
        lane_outputs+=("${lane_out}")

        if [[ ! -s "${lane_out}" ]]; then
            cutadapt \
                -j 8 \
                -a 'polyA=A{15}' \
                -e 0.1 \
                -O 10 \
                --discard-untrimmed \
                -q 20 \
                --trim-n \
                --minimum-length 25 \
                -o "${lane_out}" \
                "${raw}" \
                > "${OUT}/logs/${base}.polyA.cutadapt.log"
        fi
    done

    ((${#lane_outputs[@]} > 0)) || { echo "ERROR: no FASTQ for ${sample}" >&2; exit 1; }
    if [[ ! -s "${OUT}/polyA_reads/${sample}.polyA.fastq.gz" ]]; then
        zcat "${lane_outputs[@]}" | gzip -c > "${OUT}/polyA_reads/${sample}.polyA.fastq.gz"
    fi
done < "${OUT}/sample_ids.txt"

while read -r sample; do
    reads="${OUT}/polyA_reads/${sample}.polyA.fastq.gz"
    prefix="${OUT}/star_align/${sample}.polyA."
    bam="${prefix}Aligned.sortedByCoord.out.bam"
    ubam="${OUT}/star_align/${sample}.polyA.unique_primary.bam"

    if [[ ! -s "${bam}" ]]; then
        "${STAR_BIN}" \
            --runThreadN "${THREADS}" \
            --genomeDir "${STARREF}" \
            --readFilesIn "${reads}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${prefix}" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --quantMode GeneCounts \
            --limitBAMsortRAM 60000000000
    fi

    [[ -s "${bam}.bai" ]] || samtools index "${bam}"
    if [[ ! -s "${ubam}" ]]; then
        # STAR unique alignments have MAPQ 255; remove unmapped, secondary,
        # supplementary records (historical mask 2308 = 4+256+2048).
        # QC-failed records are NOT excluded by this historical mask.
        samtools view -@ 8 -b -q 255 -F 2308 "${bam}" > "${ubam}"
        samtools index "${ubam}"
    fi

    bedtools genomecov -ibam "${ubam}" -3 -strand + -bg \
        > "${POLY}/${sample}.polyA.plus.3end.bedgraph"
    bedtools genomecov -ibam "${ubam}" -3 -strand - -bg \
        > "${POLY}/${sample}.polyA.minus.3end.bedgraph"
done < "${OUT}/sample_ids.txt"

# Expand bedGraph intervals into per-position, per-sample observations.
for label in plus minus; do
    [[ "${label}" == plus ]] && sign='+' || sign='-'
    while read -r sample; do
        awk -v s="${sample}" -v strand="${sign}" 'BEGIN{OFS="\t"}
        {for(pos=$2; pos<$3; pos++) print $1,pos,pos+1,s,$4,strand}' \
            "${POLY}/${sample}.polyA.${label}.3end.bedgraph"
    done < "${OUT}/sample_ids.txt" |
    sort -k1,1 -k2,2n > "${POLY}/all_samples.${label}.sites.bed"

    bedtools cluster -d 25 -i "${POLY}/all_samples.${label}.sites.bed" \
        > "${POLY}/all_samples.${label}.clustered.bed"

    awk -v strand="${sign}" 'BEGIN{OFS="\t"}
    function finish_position() {
        if(!have_position) return
        n_positions++
        # Strict > preserves the earliest genomic position in an exact tie.
        if(position_sum > peak_reads){peak_reads=position_sum; peak_start=current_position}
    }
    function finish_cluster() {
        if(!have_cluster) return
        finish_position()
        print chromosome,cluster_start,cluster_end,
              "TagSeq_cluster_"strand"_"cluster_id,total_reads,strand,
              peak_start,peak_start+1,peak_reads,n_samples,n_positions
        delete observed_sample
        total_reads=0; n_samples=0; n_positions=0; peak_reads=-1; have_position=0
    }
    {
        new_cluster=$7
        if(!have_cluster || new_cluster!=cluster_id){
            finish_cluster()
            cluster_id=new_cluster; chromosome=$1; cluster_start=$2; cluster_end=$3
            have_cluster=1
        }
        if($2<cluster_start) cluster_start=$2
        if($3>cluster_end) cluster_end=$3
        if(!have_position || $2!=current_position){
            finish_position(); current_position=$2; position_sum=0; have_position=1
        }
        position_sum+=$5; total_reads+=$5
        if(!($4 in observed_sample)){observed_sample[$4]=1; n_samples++}
    }
    END{finish_cluster()}' "${POLY}/all_samples.${label}.clustered.bed" \
        > "${POLY}/all_samples.${label}.clusters.bed"
done

cat "${POLY}/all_samples.plus.clusters.bed" \
    "${POLY}/all_samples.minus.clusters.bed" |
sort -k1,1 -k2,2n > "${POLY}/all_samples.consensus_clusters.bed"

# Initial reproducibility: >=2 samples and >=10 reads across samples.
awk '$10>=2 && $5>=10' "${POLY}/all_samples.consensus_clusters.bed" \
    > "${POLY}/all_samples.initial_reproducible.bed"

# Extract the strand-aware 20-nt downstream sequence from each peak.
awk 'BEGIN{OFS="\t"}
{
    if($6=="+") print $1,$8,$8+20,$4,$5,$6
    else {s=$7-20; if(s<0)s=0; print $1,s,$7,$4,$5,$6}
}' "${POLY}/all_samples.initial_reproducible.bed" \
    > "${POLY}/all_samples.initial.downstream20.bed"

bedtools getfasta -fi "${GENOME}" \
    -bed "${POLY}/all_samples.initial.downstream20.bed" -s -name -tab \
    > "${POLY}/all_samples.initial.downstream20.tsv"

awk 'BEGIN{OFS="\t"}
{
    id=$1; seq=toupper($2); tmp=seq; a=gsub(/A/,"",tmp)
    pct=100*a/length(seq); run=(seq~/AAAAAA/)?1:0
    internal=(a>=13 || run==1)?1:0
    print id,seq,a,pct,run,internal
}' "${POLY}/all_samples.initial.downstream20.tsv" \
    > "${POLY}/all_samples.initial.internal_priming_qc.tsv"

awk -F'\t' 'FILENAME==ARGV[1]{if($6==1){split($1,x,"::"); bad[x[1]]=1}; next}!($4 in bad)' \
    "${POLY}/all_samples.initial.internal_priming_qc.tsv" \
    "${POLY}/all_samples.initial_reproducible.bed" \
    > "${POLY}/all_samples.initial.no_internal_priming.bed"

# Build 3-prime endpoints from the annotation. Column 4 is gene|transcript.
TXENDS="${POLY}/final_annotation.transcript_3ends.bed"
awk -F'\t' 'BEGIN{OFS="\t"}
$3=="transcript"{
    gid=tid=""
    if(match($9,/gene_id "[^"]+"/)){gid=substr($9,RSTART,RLENGTH); sub(/^gene_id "/,"",gid); gsub(/"/,"",gid)}
    if(match($9,/transcript_id "[^"]+"/)){tid=substr($9,RSTART,RLENGTH); sub(/^transcript_id "/,"",tid); gsub(/"/,"",tid)}
    if(gid=="" || tid=="") next
    if($7=="+") print $1,$5-1,$5,gid"|"tid,0,$7
    else print $1,$4-1,$4,gid"|"tid,0,$7
}' "${INPUT_GTF}" | sort -k1,1 -k2,2n > "${TXENDS}"

PASS="${POLY}/all_samples.initial.no_internal_priming.bed"
awk 'BEGIN{OFS="\t"}{print $1,$7,$8,$4,$5,$6}' "${PASS}" |
sort -k1,1 -k2,2n > "${POLY}/all_samples.pass.peak.bed"

# -t first makes an exact-distance tie deterministic using B-file order.
bedtools closest -a "${POLY}/all_samples.pass.peak.bed" -b "${TXENDS}" \
    -s -d -t first > "${POLY}/all_samples.pass.closest_txend.tsv"

awk 'BEGIN{OFS="\t"}{print $4,$5,$9,$10}' "${PASS}" \
    > "${POLY}/all_samples.pass.support.tsv"

awk 'BEGIN{OFS="\t"}
FILENAME==ARGV[1]{total[$1]=$2; peak[$1]=$3; samples[$1]=$4; next}
{
    if($6=="+") ext=$2-$8; else ext=$8-$2
    print $4,$10,$1,$6,$2,$3,total[$4],peak[$4],samples[$4],ext,$13
}' "${POLY}/all_samples.pass.support.tsv" \
    "${POLY}/all_samples.pass.closest_txend.tsv" \
    > "${POLY}/all_samples.pass.txend_distance.tsv"

# Conservative final support and distance filters, followed by one peak per tx.
awk -F'\t' '$2!="." && $2!="" && $9>=3 && $7>=20 && $8>=5 && $10>100 && $10<=5000' \
    "${POLY}/all_samples.pass.txend_distance.tsv" |
sort -t$'\t' -k2,2 -k9,9nr -k8,8nr -k7,7nr -k10,10n |
awk -F'\t' '!seen[$2]++' > "${POLY}/all_samples.extension.best_per_transcript.tsv"

# Gene spans used for conservative collision removal. This historical step was
# intentionally any-strand: no -s was passed to bedtools intersect.
GENES="${POLY}/final_annotation.gene_spans.bed"
awk -F'\t' 'BEGIN{OFS="\t"}
$3=="transcript"{
    gid=""
    if(match($9,/gene_id "[^"]+"/)){gid=substr($9,RSTART,RLENGTH); sub(/^gene_id "/,"",gid); gsub(/"/,"",gid)}
    if(gid=="") next
    if(!(gid in seen)){seen[gid]=1; chr[gid]=$1; start[gid]=$4-1; end[gid]=$5; strand[gid]=$7}
    if($4-1<start[gid])start[gid]=$4-1; if($5>end[gid])end[gid]=$5
}
END{for(gid in seen)print chr[gid],start[gid],end[gid],gid,0,strand[gid]}' \
    "${INPUT_GTF}" | sort -k1,1 -k2,2n > "${GENES}"

EXTBED="${POLY}/all_samples.extension.best_per_transcript.intervals.bed"
awk 'BEGIN{OFS="\t"}
FILENAME==ARGV[1]{txstart[$4]=$2; txend[$4]=$3; next}
{
    key=$2; split(key,x,"|"); gene=x[1]; transcript=x[2]
    if($4=="+"){start=txend[key]; end=$6}else{start=$5; end=txstart[key]}
    if(end>start)print $3,start,end,$1,$7,$4,gene,transcript,$9,$8,$10,key
}' "${TXENDS}" "${POLY}/all_samples.extension.best_per_transcript.tsv" |
sort -k1,1 -k2,2n > "${EXTBED}"

bedtools intersect -a "${EXTBED}" -b "${GENES}" -wa -wb \
    > "${POLY}/all_samples.extension.gene_overlaps.tsv"
awk -F'\t' '$7!=$16{print $4}' "${POLY}/all_samples.extension.gene_overlaps.tsv" |
sort -u > "${POLY}/all_samples.extension.collision_cluster.ids"

awk -F'\t' 'FILENAME==ARGV[1]{bad[$1]=1; next}!($4 in bad)' \
    "${POLY}/all_samples.extension.collision_cluster.ids" "${EXTBED}" \
    > "${POLY}/all_samples.extension.no_collision.bed"

# Remove the single historical bedtools-closest no-match placeholder (.|).
VALID="${POLY}/all_samples.extension.no_collision.valid_transcripts.bed"
awk -F'\t' '$7!="." && $8!="." && $7!="" && $8!=""' \
    "${POLY}/all_samples.extension.no_collision.bed" > "${VALID}"

# Update transcript and terminal exon only. Requires GNU awk (ARGIND).
awk -F'\t' 'BEGIN{OFS="\t"}
ARGIND==1{old_boundary[$4]=$3; next}
ARGIND==2{
    key=$7"|"$8; total_reads[key]=$5; samples[key]=$9
    if($6=="+")new_boundary[key]=$3; else new_boundary[key]=$2+1
    next
}
ARGIND==3{
    if(/^#/){print; next}
    tid=gid=""
    if(match($9,/transcript_id "[^"]+"/)){tid=substr($9,RSTART,RLENGTH); sub(/^transcript_id "/,"",tid); gsub(/"/,"",tid)}
    if(match($9,/gene_id "[^"]+"/)){gid=substr($9,RSTART,RLENGTH); sub(/^gene_id "/,"",gid); gsub(/"/,"",gid)}
    key=gid"|"tid
    if(key in new_boundary){
        if($3=="transcript"){
            if($7=="+")$5=new_boundary[key]; else $4=new_boundary[key]
            $9=$9 " tagseq_3prime \"supported\"; tagseq_samples \""samples[key]"\"; tagseq_reads \""total_reads[key]"\";"
        }else if($3=="exon"){
            if($7=="+" && $5==old_boundary[key])$5=new_boundary[key]
            else if($7=="-" && $4==old_boundary[key])$4=new_boundary[key]
        }
    }
    print
}' "${TXENDS}" "${VALID}" "${INPUT_GTF}" > "${FINAL_GTF}"

# Required integrity checks: feature counts stable and CDS byte-identical.
for file in "${INPUT_GTF}" "${FINAL_GTF}"; do
    echo "== ${file} =="
    awk -F'\t' '$3=="transcript"{t++}$3=="exon"{e++}$3=="CDS"{c++}
        END{print "transcripts",t+0; print "exons",e+0; print "CDS",c+0}' "${file}"
done

input_cds_md5=$(awk -F'\t' '$3=="CDS"' "${INPUT_GTF}" | md5sum | awk '{print $1}')
final_cds_md5=$(awk -F'\t' '$3=="CDS"' "${FINAL_GTF}" | md5sum | awk '{print $1}')
[[ "${input_cds_md5}" == "${final_cds_md5}" ]] || {
    echo "ERROR: CDS records changed during TagSeq refinement" >&2
    exit 1
}

awk -F'\t' '$3=="transcript" && $9~/tagseq_3prime "supported"/{n++}
END{print "TagSeq-modified transcripts:",n+0}' "${FINAL_GTF}"

if command -v cellranger >/dev/null; then
    cellranger mkgtf "${FINAL_GTF}" "${CHECKED_GTF}"
else
    echo "NOTE: cellranger not in PATH; skipped mkgtf validation." >&2
fi

echo "Final TagSeq-refined GTF: ${FINAL_GTF}"
