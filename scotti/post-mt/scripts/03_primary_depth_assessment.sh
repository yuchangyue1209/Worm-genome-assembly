#!/usr/bin/env bash
# Environment: scotti_purge_eval; minimap2 absolute path avoids environment ambiguity.
# Assessment only: DOES NOT purge sequences or choose cutoffs automatically.
source "$(dirname "$0")/00_config.sh"
OUT="${OUT:-$ROOT/qc/primary_purge_assessment_reconstructed}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"; cd "$OUT"
/home/cyu/.conda/envs/genome_assembly/bin/minimap2 -x map-hifi -t 20 "$PRIMARY" "$HIFI" \
  2> minimap2.log | gzip -c > primary_vs_hifi.paf.gz
pbcstat primary_vs_hifi.paf.gz 2>&1 | tee pbcstat.log
# PB.stat: coverage-bin distribution; inspect together with BUSCO and self-alignment.
