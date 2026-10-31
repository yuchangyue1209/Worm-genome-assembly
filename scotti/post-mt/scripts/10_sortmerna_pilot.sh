#!/usr/bin/env bash
# Environment: rrna_filter. Reuses previously prepared first-100,000-pair pilot.
source "$(dirname "$0")/00_config.sh"
QC="$ANNO/results/rnaseq_qc"
REF="$QC/02_overrepresented_check/solidus_nuclear_rrna.fa"
OUT="${OUT:-$QC/04_sortmerna_pilot_reconstructed}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"
sortmerna --version > "$OUT/version.txt" 2>&1
sortmerna --ref "$REF" \
  --reads "$QC/03_rrna_pilot/pilot_R1.fastq" \
  --reads "$QC/03_rrna_pilot/pilot_R2.fastq" \
  --workdir "$OUT/work" --aligned "$OUT/rrna" --other "$OUT/non_rrna" \
  --fastx --paired_in --out2 --threads 12 2>&1 | tee "$OUT/run.log"
# paired_in: if either mate hits, both go to the rRNA-associated output.
