#!/usr/bin/env bash
# Environment: genome_assembly. Historical lineage: metazoa_odb10, 954 groups.
source "$(dirname "$0")/00_config.sh"
OUT="${OUT:-$ROOT/qc/primary_busco_reconstructed}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"
busco -i "$PRIMARY" -o primary_metazoa_odb10 --out_path "$OUT" \
  -l metazoa_odb10 -m genome --cpu 20 2>&1 | tee "$OUT/run.log"
