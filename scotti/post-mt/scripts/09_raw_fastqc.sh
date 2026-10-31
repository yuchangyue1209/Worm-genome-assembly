#!/usr/bin/env bash
# Current rnaseq_annotation executable uses --outdir; -o was rejected.
source "$(dirname "$0")/00_config.sh"
OUT="${OUT:-$ANNO/results/rnaseq_qc/01_raw_fastqc_reconstructed}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"
fastqc -t 2 --outdir "$OUT" "$R1" "$R2" 2>&1 | tee "$OUT/run.log"
