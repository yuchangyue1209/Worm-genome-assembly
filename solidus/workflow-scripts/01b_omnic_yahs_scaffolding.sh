#!/usr/bin/env bash
# Reconstructed runnable workflow; not a recovered historical command log.
# Usage: bash 01b_omnic_yahs_scaffolding.sh CONTIGS.fa R1.fq.gz R2.fq.gz NEW_OUTDIR
# Or:    bash 01b_omnic_yahs_scaffolding.sh --bam CONTIGS.fa ALIGNMENTS.bam NEW_OUTDIR
# --bam requires alignments against these exact contigs. Duplicates must already
# be removed or marked; the script filters marked duplicates and name-sorts.
set -euo pipefail

MODE=reads
if [[ "${1:-}" == --bam ]]; then
    [[ $# -eq 4 ]] || { echo 'Usage: --bam CONTIGS.fa ALIGNMENTS.bam NEW_OUTDIR' >&2; exit 2; }
    MODE=bam; ASM=$2; INPUT_BAM=$3; OUT=$4
else
    [[ $# -eq 4 ]] || { echo 'Usage: CONTIGS.fa R1.fq.gz R2.fq.gz NEW_OUTDIR' >&2; exit 2; }
    ASM=$1; R1=$2; R2=$3; OUT=$4
fi
THREADS=${THREADS:-16}
SORT_THREADS=${SORT_THREADS:-4}
SORT_MEM=${SORT_MEM:-1G} # per sorting thread
for tool in samtools yahs; do
    command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
[[ -s "$ASM" ]] || { echo "Missing contig FASTA: $ASM" >&2; exit 1; }
if [[ $MODE == reads ]]; then
    command -v bwa-mem2 >/dev/null || { echo 'Missing bwa-mem2' >&2; exit 1; }
    for f in "$R1" "$R2"; do [[ -s "$f" ]] || { echo "Missing reads: $f" >&2; exit 1; }; done
else
    [[ -s "$INPUT_BAM" ]] || { echo "Missing BAM: $INPUT_BAM" >&2; exit 1; }
    samtools quickcheck -v "$INPUT_BAM"
fi
[[ ! -e "$OUT" ]] || { echo "Output already exists; use a new directory: $OUT" >&2; exit 1; }
mkdir -p "$OUT"
cp "$ASM" "$OUT/contigs.fa"
samtools faidx "$OUT/contigs.fa"
{
    echo 'Reconstructed workflow: q10, resolution 500000, no contig error correction.'
    echo 'Mapping and duplicate-processing choices are not verified historical settings.'
    samtools --version
    yahs --version
    if [[ $MODE == reads ]]; then bwa-mem2 version; fi
} > "$OUT/software_versions.txt" 2>&1

if [[ $MODE == reads ]]; then
    bwa-mem2 index "$OUT/contigs.fa" > "$OUT/index.log" 2>&1
    # Treat Omni-C mates independently; preserve the 5-prime alignment for
    # split reads. These are reconstructed mapping choices, not recovered flags.
    bwa-mem2 mem -5SP -t "$THREADS" "$OUT/contigs.fa" "$R1" "$R2" \
        2> "$OUT/mapping.log" |
        samtools view -u -F 2308 - |
        samtools sort -n -@ "$SORT_THREADS" -m "$SORT_MEM" \
            -T "$OUT/name_tmp" -o "$OUT/primary.namesort.bam" -
    samtools fixmate -m -@ "$SORT_THREADS" "$OUT/primary.namesort.bam" "$OUT/fixmate.bam"
    samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "$OUT/coord_tmp" \
        -o "$OUT/fixmate.coordsort.bam" "$OUT/fixmate.bam"
    samtools markdup -r -s -@ "$SORT_THREADS" \
        "$OUT/fixmate.coordsort.bam" "$OUT/deduplicated.bam" \
        2> "$OUT/markdup.log"
    INPUT_BAM=$OUT/deduplicated.bam
fi
# Filter unmapped, secondary, supplementary and marked duplicate records.
# Name sorting ensures that YaHS applies its -q mapping-quality threshold.
samtools view -u -F 3332 "$INPUT_BAM" |
    samtools sort -n -@ "$SORT_THREADS" -m "$SORT_MEM" \
        -T "$OUT/yahs_name_tmp" -o "$OUT/omnic.yahs.namesort.bam" -
samtools quickcheck -v "$OUT/omnic.yahs.namesort.bam"
samtools flagstat "$OUT/omnic.yahs.namesort.bam" > "$OUT/alignments.flagstat.txt"
yahs -q 10 -r 500000 --no-contig-ec -o "$OUT/ssolidus_q10_r500k_no_contig_ec" \
    "$OUT/contigs.fa" "$OUT/omnic.yahs.namesort.bam" > "$OUT/yahs.log" 2>&1
PREFIX=$OUT/ssolidus_q10_r500k_no_contig_ec
[[ -s "${PREFIX}_scaffolds_final.fa" && -s "${PREFIX}_scaffolds_final.agp" ]] || {
    echo "Missing final scaffolds; inspect $OUT/yahs.log" >&2; exit 1;
}
samtools faidx "${PREFIX}_scaffolds_final.fa"
printf 'Scaffolds: %s\nAGP: %s\n' "${PREFIX}_scaffolds_final.fa" "${PREFIX}_scaffolds_final.agp"
# Contact-map review, mitochondrial separation and FCS-GX curation are later
# stages. This script does not overwrite or recreate the curated final release.
