#!/usr/bin/env bash
# Build a STAR index, align the 15 PRJNA304161 paired libraries and generate
# AUGUSTUS intron hints.
set -euo pipefail

ROOT=/work/cyu/annotation/ssol-annotation
GENOME=/work/cyu/assembly/worm/hifiasm/purge_q10/worm.q10_r500k_no_contig_ec.nuclear.final.fa
RNASEQ=/work/shared/transcriptomes/hebert-2016/PRJNA304161/SRP066813
THREADS=32
INDEX=${ROOT}/results/rnaseq/star_index
ALIGN=${ROOT}/results/rnaseq/star_align
TMP=${ROOT}/temp/star
BAM=${ALIGN}/Aligned.sortedByCoord.out.bam
mkdir -p "${INDEX}" "${ALIGN}" "${TMP}/sort"

for cmd in STAR samtools bam2hints; do
    command -v "${cmd}" >/dev/null || { echo "ERROR: ${cmd} is not in PATH" >&2; exit 1; }
done

shopt -s nullglob
r1=("${RNASEQ}"/*_1.fastq.gz)
r2=("${RNASEQ}"/*_2.fastq.gz)
shopt -u nullglob
(( ${#r1[@]} > 0 && ${#r1[@]} == ${#r2[@]} )) || {
    echo "ERROR: paired FASTQ files are missing or mismatched" >&2; exit 1;
}
mapfile -t r1 < <(printf '%s\n' "${r1[@]}" | sort)
mapfile -t r2 < <(printf '%s\n' "${r2[@]}" | sort)
r1_csv=$(IFS=,; echo "${r1[*]}")
r2_csv=$(IFS=,; echo "${r2[*]}")

genome_size=$(awk '/^[^>]/{n+=length($0)} END{print n}' "${GENOME}")
sa_bases=$(python3 -c "import math; print(min(14, int(math.log2(${genome_size})/2 - 1)))")

if [[ ! -s "${INDEX}/Genome" ]]; then
    STAR --runMode genomeGenerate --runThreadN "${THREADS}" \
        --genomeDir "${INDEX}" --genomeFastaFiles "${GENOME}" \
        --genomeSAindexNbases "${sa_bases}"
fi

STAR --runMode alignReads --runThreadN "${THREADS}" --genomeDir "${INDEX}" \
    --readFilesIn "${r1_csv}" "${r2_csv}" --readFilesCommand zcat \
    --twopassMode Basic --sjdbOverhang 149 \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 20 --outFileNamePrefix "${ALIGN}/"

samtools sort -@ "${THREADS}" -T "${TMP}/sort/st" \
    -o "${BAM}" "${ALIGN}/Aligned.out.bam"
samtools index "${BAM}"
rm -f "${ALIGN}/Aligned.out.bam"

# GeneMark-ETP calls StringTie, which needs XS strand tags on spliced STAR
# alignments. Abort here rather than allowing a later parse_set.pl division-by-zero.
if ! samtools view "${BAM}" |
    awk '$6 ~ /N/ {for(i=12;i<=NF;i++) if($i ~ /^XS:A:/) found=1; if(found) exit} END{exit(found ? 0 : 1)}'
then
    echo "ERROR: no XS:A: tag was found on a spliced alignment" >&2
    exit 1
fi

bam2hints --intronsonly --in="${BAM}" --out="${ROOT}/results/rnaseq/hints.gff"
cat "${ALIGN}/Log.final.out"
wc -l "${ROOT}/results/rnaseq/hints.gff"
