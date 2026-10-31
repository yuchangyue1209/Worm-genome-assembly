#!/usr/bin/env bash
# Environment: genome_assembly. BLAST evidence; no automatic removal of local mt hits.
source "$(dirname "$0")/00_config.sh"
OUT="${OUT:-$ROOT/qc/mt_source_tracking_reconstructed}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"
for label in primary ragtag samba; do
  case "$label" in primary) ASM="$PRIMARY";; ragtag) ASM="$RAGTAG";; samba) ASM="$NEW_SAMBA";; esac
  blastn -task blastn -query "$MT" -subject "$ASM" -evalue 1e-10 \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
    -out "$OUT/$label.mt_hits.tsv"
done
# Review union of aligned intervals, not summed HSP lengths. Keep large nuclear sequences
# containing local mt-like segments; those are not automatically mt-only contigs.
