#!/usr/bin/env bash
# Environment: genome_assembly. Select one documented evaluation; never overwrites results.
# Usage: bash 20_busco_proteins.sh raw10|longest10|longest12
source "$(dirname "$0")/00_config.sh"
case "${1:-}" in
 raw10) PROT="$ANNO/results/braker3_etp_v1_retry2/braker.aa"; DIR=01_raw_braker_busco; NAME=raw_braker_metazoa; LINEAGE=metazoa_odb10;;
 longest10) PROT="$ANNO/results/annotation_qc/02_longest_isoform_busco/Scotti_braker.longest.proteins.fa"; DIR=02_longest_isoform_busco; NAME=longest_metazoa; LINEAGE=metazoa_odb10;;
 longest12) PROT="$ANNO/results/annotation_qc/02_longest_isoform_busco/Scotti_braker.longest.proteins.fa"; DIR=05_longest_busco_odb12; NAME=scotti_longest_metazoa_odb12; LINEAGE=metazoa_odb12;;
 *) echo 'Usage: raw10|longest10|longest12'; exit 1;;
esac
OUT="$ANNO/results/annotation_qc/$DIR"
[[ -s "$PROT" && ! -e "$OUT/$NAME" ]]
mkdir -p "$OUT"
busco -i "$PROT" -m proteins -l "$LINEAGE" --cpu 20 --out_path "$OUT" -o "$NAME" 2>&1 | tee "$OUT/run.log"
