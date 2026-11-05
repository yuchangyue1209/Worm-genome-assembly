#!/usr/bin/env bash
# Environment: genome_assembly. Reuses the already extracted 72 regions.
# CRITICAL: -c is required before interpreting PAF identity quantitatively.
source "$(dirname "$0")/00_config.sh"
OUT="$ANNO/results/annotation_qc/04_busco_candidate_mapping"
[[ -s "$OUT/primary_candidate_regions.fa" ]]
PAF="$OUT/candidates_vs_nuclear.baselevel.paf"
[[ ! -e "$PAF" ]] || { echo 'Base-level result exists; do not overwrite'; exit 1; }
minimap2 -x asm5 -c --cs=short --secondary=yes -N 10 -t 16 \
 "$NUCLEAR" "$OUT/primary_candidate_regions.fa" > "$PAF" 2> "$OUT/minimap2.baselevel.log"
# Earlier no-CIGAR PAF matches/lengths were approximate; its 37/35 classification was withdrawn.
