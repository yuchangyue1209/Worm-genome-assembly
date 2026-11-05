#!/usr/bin/env bash
# Environment: rnaseq_annotation. Command supplied; completion not reported.
source "$(dirname "$0")/00_config.sh"
READS="$ANNO/data/rnaseq/rrna_screened_v1/non_rrna"
OUT="$ANNO/results/rnaseq_qc/06_filtered_fastqc_chunk00001"
mkdir "$OUT"
fastqc -t 2 --outdir "$OUT" "$READS/chunk00001_R1.fastq.gz" \
  "$READS/chunk00001_R2.fastq.gz" 2>&1 | tee "$OUT/run.log"
# This first-chunk check does not establish whole-library post-filter quality.
