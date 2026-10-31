#!/usr/bin/env bash
# Environment: genome_assembly with blastn, makeblastdb and Python.
# Reconstruct the initial mt screen. Use a new output directory, not the old qc/mt results.
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
need "$RAW_ASM" "$SOL_REF"
OUT="$ROOT/qc/mt_reproduced"
fresh "$OUT"
makeblastdb -in "$RAW_ASM" -dbtype nucl -out "$OUT/scotti_db"
blastn -task blastn -query "$SOL_REF" -db "$OUT/scotti_db" -evalue 1e-10 -max_target_seqs 100000 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovhsp' -out "$OUT/solidus_mt_vs_scotti.tsv"
sort -k4,4nr "$OUT/solidus_mt_vs_scotti.tsv" > "$OUT/solidus_mt_vs_scotti.sorted.tsv"
python "$HELPER" summarize "$OUT/solidus_mt_vs_scotti.tsv" > "$OUT/mt_hits_summary.tsv"
# Four IDs selected after manual review; this is NOT a general similarity-only deletion rule.
# Writes a reproduced provisional split, not a final nuclear release.
python "$HELPER" split-candidates "$RAW_ASM" "$OUT"
blastn -task blastn -query "$OUT/scotti_mt_candidates.fa" -subject "$OUT/scotti_mt_candidates.fa" -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore' -out "$OUT/mt_candidates_self.tsv"
awk '$1!=$2 || $5!=$7 || $6!=$8' "$OUT/mt_candidates_self.tsv" | sort -k11,11nr > "$OUT/mt_candidates_self.filtered.tsv"
# Existing extraction input remains $CAND, made earlier from the same four candidates.
