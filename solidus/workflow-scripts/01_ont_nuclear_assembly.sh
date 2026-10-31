#!/usr/bin/env bash
# Historical ONT nuclear-assembly command record. Run the documented stages
# individually after checking paths; this file preserves exploratory branches.

# ONT assembly

filtlong --min_length 1000 --keep_percent 95 worm_all.fastq.gz \
  | gzip > worm_all.filt1k.top95.fastq.gz


flye --nano-raw worm_all.filt1k.top95.fastq.gz \
     --genome-size 550m \
     --out-dir flye_top95 \
     --threads 32


test
# Assembly contiguity

test env
quast -t 16 -o quast.flye flye_top95/assembly.fasta

# Assembly completeness with the metazoa or eukaryota BUSCO lineage
busco -i flye_top95/assembly.fasta -l metazoa_odb10 -m genome -c 32 -o busco.flye



#95% qc hifiasm ont
cd /work/cyu/hifiasm

./hifiasm --ont \
  -o /work/cyu/assembly/worm/hifiasm/worm.asm \
  -t32 \
  /work/cyu/assembly/worm/worm_all.filt1k.top95.fastq.gz \
  2> /work/cyu/assembly/worm/hifiasm/worm.asm.log

cd /work/cyu/assembly/worm/hifiasm

# QUAST
quast -t 32 -o quast.worm_hifiasm worm.asm.p_ctg.fa

# BUSCO using the metazoa_odb10 lineage
busco -i worm.asm.p_ctg.fa \
      -l metazoa_odb10 \
      -m genome \
      -c 32 \
      -o busco.worm_hifiasm



cd /work/cyu/assembly/worm/hifiasm

# 1. Align ONT reads to the primary assembly
minimap2 -x map-ont -t 32 \
  worm.asm.p_ctg.fa \
  /work/cyu/assembly/worm/worm_all.filt1k.top95.fastq.gz \
  > aln.r1.paf

# 2. First round of Racon polishing
racon -t 32 \
  /work/cyu/assembly/worm/worm_all.filt1k.top95.fastq.gz \
  aln.r1.paf \
  worm.asm.p_ctg.fa \
  > worm.polish.r1.fa

















#90% min q10
filtlong --keep_percent 90 --min_mean_q 10 ../worm_all.fastq.gz \
  | gzip > worm_all.top90.q10.fastq.gz

./hifiasm --ont \
  -o /work/cyu/assembly/worm/hifiasm/worm.asm \
  -t32 \
  /work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz \
  2> /work/cyu/assembly/worm/hifiasm/worm.asm.log


quast -t 16 \
  -o /work/cyu/assembly/worm/hifiasm/quast.worm_hifiasm_q10 \
  /work/cyu/assembly/worm/hifiasm/worm_q10.asm.p_ctg.fa


busco -i worm.asm.p_ctg.fa \
      -l metazoa_odb10 \
      -m genome \
      -c 32 \
      -o busco.worm_hifiasm


cd /work/cyu/assembly/worm/hifiasm

# Round 1: read alignment
minimap2 -x map-ont -t 32 \
  worm_q10.asm.p_ctg.fa \
  /work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz \
  > worm_q10.r1.paf

# Round 1: Racon
racon -t 32 \
  /work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz \
  worm_q10.r1.paf \
  worm_q10.asm.p_ctg.fa \
  > worm_q10.polish.r1.fa

cd /work/cyu/assembly/worm/hifiasm

# Realign reads to the round-1 polished assembly
minimap2 -x map-ont -t 32 \
  worm_q10.polish.r1.fa \
  /work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz \
  > worm_q10.r2.paf

# Second round of Racon polishing
racon -t 32 \
  /work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz \
  worm_q10.r2.paf \
  worm_q10.polish.r1.fa \
  > worm_q10.polish.r2.fa

#medaka deep learning based polishing medaka env
medaka_consensus \
  -i /work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz \
  -d worm_q10.polish.r2.fa \
  -o medaka_q10 \
  -t 32 \
  -m r1041_e82_400bps_sup_v5.2.0

busco \
  -i worm_q10.polish.medaka.fa \
  -l metazoa_odb10 \
  -o busco.worm_medaka \
  -m genome \
  -c 32
quast -t 16 \
  -o quast.worm_medaka \
  worm_q10.polish.medaka.fa
cd /work/cyu/assembly/worm/hifiasm

# Define convenient links for the assembly and filtered ONT reads
ln -s worm_q10.polish.medaka.fa asm.fa
ln -s /work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz ont.q10.fastq.gz


#run purge
conda activate test

cd /work/cyu/assembly/worm/hifiasm
mkdir -p purge_q10
cd purge_q10

# Path variables used below
ASM=../worm_q10.polish.medaka.fa
READS=/work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz
PD=/work/cyu/purge_dups/src
# Align ONT reads to the polished assembly with minimap2

minimap2 -t 32 -x map-ont "$ASM" "$READS" | gzip -c > reads_vs_asm.paf.gz


conda activate poolseq_env
# Calculate the coverage distribution
"$PD"/pbcstat reads_vs_asm.paf.gz
# Expected outputs: PB.base.cov and PB.stat
ls PB.*
# Estimate purge_dups cutoffs from PB.stat
"$PD"/calcuts PB.stat > cutoffs
cat cutoffs

# The following commands require the poolseq_env environment with purge_dups
conda activate poolseq_env

# 2.1 Calculate the coverage distribution
"$PD"/pbcstat reads_vs_asm.paf.gz

# Expected outputs:
#   PB.base.cov
#   PB.stat
ls PB.*

# 2.2 Estimate cutoffs from PB.stat
"$PD"/calcuts PB.stat > cutoffs

echo "==== cutoffs ===="
cat cutoffs
echo "================="

# Ensure that purge_dups is available in the active environment

ASM=/work/cyu/assembly/worm/hifiasm/worm_q10.polish.medaka.fa
PD=/work/cyu/purge_dups/src

"$PD"/split_fa "$ASM" > asm.split.fa
ls asm.split.fa

# Assembly self-alignment requires an environment containing minimap2
conda activate test  # Alternatively use another environment with minimap2

minimap2 -t 32 -xasm5 -DP asm.split.fa asm.split.fa \
  | gzip -c > asm.split.self.paf.gz

ls asm.split.self.paf.gz


conda activate poolseq_env  # Return to the environment containing purge_dups

"$PD"/purge_dups -2 \
  -T cutoffs \
  -c PB.base.cov \
  asm.split.self.paf.gz > dups.bed

head dups.bed

"$PD"/get_seqs dups.bed "$ASM"

ls purged.fa hap.fa



mv purged.fa worm_q10.medaka.purged.fa
mv hap.fa    worm_q10.medaka.hap.fa

ls -lh worm_q10.medaka*.fa






# Test the asm10 preset using a more permissive self-alignment
conda activate test
minimap2 -t 32 -xasm10 -DP asm.split.fa asm.split.fa \
  | gzip -c > asm.split.self.asm10.paf.gz

conda activate poolseq_env
$PD/purge_dups -2 \
  -T cutoffs \
  -c PB.base.cov \
  asm.split.self.asm10.paf.gz > dups.asm10.bed

$PD/get_seqs dups.asm10.bed "$ASM"

mv purged.fa worm_q10.medaka.purged.asm10.fa
mv hap.fa    worm_q10.medaka.hap.asm10.fa

# Run QUAST and BUSCO on the new purged assembly
quast -t 16 -o quast.worm_medaka_purged.asm10 worm_q10.medaka.purged.asm10.fa
busco -i worm_q10.medaka.purged.asm10.fa -l metazoa_odb10 -m genome -c 32 \
      -o busco.worm_medaka_purged.asm10



# The result was identical to the asm5 run.




minimap2 -t 32 -ax map-ont --cs -p 0.3 \
  "$ASM" "$READS" \
  | samtools sort -@ 16 -o ont_vs_asm.chimera.bam

samtools index ont_vs_asm.chimera.bam

samtools view ont_vs_asm.chimera.bam \
  | awk '
    $0 ~ /SA:Z:/ {
      chr=$3;
      split_count[chr]++
    }
    END {
      for (c in split_count) {
        print c, split_count[c]
      }
  }' | sort -k2,2nr | head





#!/usr/bin/env bash
set -euo pipefail

# Map ONT reads back to assembly for chimera/coverage-based QC
# Outputs sorted BAM + index

T="${T:-32}"

ASM="${ASM:-worm_q10.medaka.purged.asm10.fa}"
READS="${READS:-/work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz}"
OUTBAM="${OUTBAM:-ont_vs_asm.chimera.bam}"

module_loaded=false

echo "[01] Mapping ONT reads -> assembly"
minimap2 -t "${T}" -ax map-ont --cs -p 0.3 "${ASM}" "${READS}" \
  | samtools sort -@ 16 -o "${OUTBAM}"

samtools index "${OUTBAM}"
echo "[01] Done: ${OUTBAM} (+ .bai)"


#!/usr/bin/env bash
set -euo pipefail

# Compute per-contig mean depth (mosdepth regions), length, and GC.
# Identify suspect short contigs with depth outliers.

T="${T:-16}"

ASM="${ASM:-worm_q10.medaka.purged.asm10.fa}"
BAM="${BAM:-ont_vs_asm.chimera.bam}"

PREFIX="${PREFIX:-asm_cov}"              # mosdepth prefix
OUT_LEN="${OUT_LEN:-contig.len.tsv}"
OUT_GC="${OUT_GC:-contig.gc.tsv}"
OUT_DEPTH="${OUT_DEPTH:-contig.len.depth.tsv}"
OUT_ALL="${OUT_ALL:-contig.len.depth.gc.tsv}"

# Outlier thresholds (tune if needed)
MIN_LEN_SUSPECT="${MIN_LEN_SUSPECT:-0}"
MAX_LEN_SUSPECT="${MAX_LEN_SUSPECT:-50000}"     # 50 kb
LOW_DEPTH_MULT="${LOW_DEPTH_MULT:-0.25}"





#
cd /work/cyu/assembly/worm/hifiasm/purge_q10

ASM=worm_q10.medaka.purged.asm10.fa
BAM=ont_vs_asm.chimera.bam
T=32

# 1. Ensure that the FASTA index exists
samtools faidx "$ASM"

# 2. Extract contig lengths
cut -f1,2 "$ASM.fai" > contig.len.tsv

# 3. Calculate mean contig depth from the existing BAM
samtools coverage -w 0 "$BAM" > contig.coverage.tsv

# 4) contig GC%
seqkit fx2tab -n -g "$ASM" > contig.gc.tsv

# 5. Combine contig length, mean depth, and GC into one table
awk 'BEGIN{OFS="\t"}
  NR==FNR{len[$1]=$2; next}
  FNR==1{next}
  {depth[$1]=$7}  # Column 7 is mean depth in this samtools coverage output
  END{
    for(c in len) print c, len[c], (c in depth?depth[c]:"NA")
  }' contig.len.tsv contig.coverage.tsv \
  | sort -k2,2nr > contig.len.depth.tsv

# Add GC content
awk 'BEGIN{OFS="\t"}
  NR==FNR{gc[$1]=$2; next}
  {print $0, ( $1 in gc ? gc[$1] : "NA")}' contig.gc.tsv contig.len.depth.tsv \
  > contig.len.depth.gc.tsv

head contig.len.depth.gc.tsv



awk 'BEGIN{OFS="\t"} $3!="NA"{print $3}' contig.len.depth.gc.tsv \
  | sort -n \
  | awk '{
      a[NR]=$1
    }
    END{
      if(NR==0){print "NA"; exit}
      if(NR%2==1) print a[(NR+1)/2];
      else print (a[NR/2]+a[NR/2+1])/2
    }'



MED=50   # Replace with the median depth calculated in the preceding step

awk -v med="$MED" 'BEGIN{OFS="\t"}
  {
    contig=$1; len=$2; depth=$3; gc=$4
    if(depth=="NA") next
    ratio=depth/med
    if(ratio<0.3 || ratio>3){
      print contig, len, depth, ratio, gc
    }
  }' contig.len.depth.gc.tsv \
  | sort -k4,4nr | head -50


cd /work/cyu/assembly/worm/hifiasm/purge_q10

ls -lh worm_q10.medaka.purged.asm10.fa



#decotamination
#!/bin/bash
set -euo pipefail

############################################
# Contamination screening and rDNA removal
############################################

cd /work/cyu/assembly/worm/hifiasm/purge_q10

ASM=worm_q10.medaka.purged.asm10.fa
READS=/work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz
THREADS=32

############################################
# 1. Map ONT reads back to assembly
############################################

minimap2 -t ${THREADS} \
    -ax map-ont \
    --cs \
    -p 0.3 \
    ${ASM} ${READS} \
| samtools sort -@16 \
    -o ont_vs_asm.chimera.bam

samtools index ont_vs_asm.chimera.bam

############################################
# 2. Calculate mean coverage for each contig
############################################

samtools faidx ${ASM}

awk 'BEGIN{OFS="\t"}{print $1,0,$2}' \
    ${ASM}.fai \
    > contigs.bed

mosdepth \
    -t 16 \
    --by contigs.bed \
    asm_cov \
    ont_vs_asm.chimera.bam

zcat asm_cov.regions.bed.gz \
| awk 'BEGIN{OFS="\t"}{print $1,$4}' \
> contig.depth.tsv

############################################
# 3. Calculate contig length and GC
############################################

cut -f1,2 ${ASM}.fai \
> contig.len.tsv

seqkit fx2tab -n -g ${ASM} \
> contig.gc.tsv

############################################
# 4. Merge coverage / GC / length
############################################

awk 'BEGIN{OFS="\t"}
NR==FNR{
    len[$1]=$2
    next
}
{
    depth[$1]=$2
}
END{
    for(c in len){
        print c,len[c],depth[c]
    }
}' \
contig.len.tsv \
contig.depth.tsv \
> tmp.len.depth.tsv

awk 'BEGIN{OFS="\t"}
NR==FNR{
    gc[$1]=$2
    next
}
{
    print $1,$2,$3,gc[$1]
}' \
contig.gc.tsv \
tmp.len.depth.tsv \
| sort -k2,2nr \
> contig.len.depth.gc.tsv

############################################
# 5. Identify short depth-outlier contigs
############################################

MED=$(awk '{print $3}' contig.len.depth.gc.tsv \
      | sort -n \
      | awk '{
          a[NR]=$1
      }
      END{
          print a[int((NR+1)/2)]
      }')

awk -v med=$MED '
BEGIN{OFS="\t"}
{
    ratio=$3/med

    if($2<50000 && (ratio<0.3 || ratio>3))
        print
}' \
contig.len.depth.gc.tsv \
> suspect.lt50k.depthOutliers.tsv

cut -f1 \
suspect.lt50k.depthOutliers.tsv \
> suspect.contigs.list

############################################
# 6. Extract suspect contigs
############################################

seqkit grep \
    -f suspect.contigs.list \
    ${ASM} \
> suspect.fa

############################################
# 7. rRNA screening
############################################

barrnap \
    --kingdom bac \
    suspect.fa \
> suspect.bac_rRNA.gff

barrnap \
    --kingdom euk \
    suspect.fa \
> suspect.euk_rRNA.gff

grep -v "^#" suspect.bac_rRNA.gff \
| cut -f1 \
| sort -u \
> suspect.bac_rRNA.contigs.list

grep -v "^#" suspect.euk_rRNA.gff \
| cut -f1 \
| sort -u \
> suspect.euk_rRNA.contigs.list

############################################
# 8. Compare bacterial/eukaryotic hits
############################################

sort suspect.bac_rRNA.contigs.list \
> bac.list

sort suspect.euk_rRNA.contigs.list \
> euk.list

comm -12 bac.list euk.list \
> both.list

comm -23 bac.list euk.list \
> bac_only.list

comm -13 bac.list euk.list \
> euk_only.list

############################################
# 9. Host contamination screen
############################################

diamond blastx \
    -d /work/cyu/diamond/stickleback.dmnd \
    -q suspect.fa \
    -o suspect.vs_stickle.m8.tsv \
    -p ${THREADS} \
    -e 1e-10 \
    -k 1 \
    --sensitive \
    --outfmt 6 \
    qseqid \
    sseqid \
    pident \
    length \
    evalue \
    bitscore \
    stitle

############################################
# 10. Remove rDNA-associated contigs
############################################

seqkit grep \
    -f suspect.euk_rRNA.contigs.list \
    ${ASM} \
> worm.rDNA_contigs.fa

seqkit grep \
    -v \
    -f suspect.euk_rRNA.contigs.list \
    ${ASM} \
> worm.no_rDNA.fa

############################################
# 11. Evaluate final assembly
############################################

quast \
    -t 16 \
    -o quast.worm.no_rDNA \
    worm.no_rDNA.fa

busco \
    -i worm.no_rDNA.fa \
    -l metazoa_odb10 \
    -m genome \
    -c 32 \
    -o busco.worm.no_rDNA

############################################
# 12. Produce ranked contig table
############################################

samtools faidx worm.no_rDNA.fa

cut -f1,2 worm.no_rDNA.fa.fai \
| sort -k2,2nr \
| awk '
BEGIN{
    OFS="\t"
    print "rank","contig","length_bp","length_Mb"
}
{
    printf "%d\t%s\t%d\t%.6f\n",
        NR,$1,$2,$2/1000000
}' \
> worm.no_rDNA.contig_lengths.with_rank.tsv
