#!/usr/bin/env bash
# Environment: genome_assembly. Historical SAMBA validation experiment.
source "$(dirname "$0")/00_config.sh"
OUT="${OUT:-$ROOT/qc/samba_hifi_read_support_reconstructed}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"
minimap2 -ax map-hifi --secondary=no -t 20 "$NEW_SAMBA" "$HIFI" \
  2> "$OUT/minimap2.log" | samtools sort -@ 8 -m 1G \
  -T "$OUT/sort_tmp" -o "$OUT/hifi_vs_new_samba.bam" -
samtools index -@ 8 "$OUT/hifi_vs_new_samba.bam"
samtools quickcheck -v "$OUT/hifi_vs_new_samba.bam"
samtools flagstat -@ 8 "$OUT/hifi_vs_new_samba.bam" | tee "$OUT/flagstat.txt"
# Whole-assembly mapping rate does not establish junction correctness.
