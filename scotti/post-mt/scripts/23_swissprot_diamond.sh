#!/usr/bin/env bash
# Reconstructed completed search on ALL 11,471 BRAKER isoforms, not just longest proteins.
source "$(dirname "$0")/00_config.sh"
PROT="$ANNO/results/braker3_etp_v1_retry2/braker.aa"
OUT="$ANNO/results/functional/01_swissprot"
DB=/work/cyu/annotation/ssol-annotation/results/functional_final/swissprot/uniprot_sprot
DBFA=/work/cyu/annotation/ssol-annotation/results/functional_final/database/uniprot_sprot.fasta
DIAMOND=$(command -v diamond || true)
if [[ -z "$DIAMOND" ]];then
 for env in genome_annotation eggnog_mapper braker3_annotation genome_assembly;do
  candidate="/home/cyu/.conda/envs/$env/bin/diamond"
  if [[ -x "$candidate" ]];then DIAMOND="$candidate";break;fi
 done
fi
[[ -n "$DIAMOND" && -s "$PROT" ]]
mkdir -p "$(dirname "$OUT")";mkdir "$OUT"
"$DIAMOND" version > "$OUT/diamond.version.txt"
sha256sum "$PROT" > "$OUT/query.sha256"
if [[ ! -s "${DB}.dmnd" ]];then
 [[ -s "$DBFA" ]];DB="$OUT/uniprot_sprot";"$DIAMOND" makedb --in "$DBFA" --db "$DB"
fi
printf '%s\n' "$DB" > "$OUT/database.path.txt"
"$DIAMOND" blastp --db "$DB" --query "$PROT" --out "$OUT/scotti_vs_swissprot.top5.tsv" \
 --threads 32 --sensitive --evalue 1e-5 --max-target-seqs 5 \
 --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle 2>&1 | tee "$OUT/run.log"
LC_ALL=C sort -t $'\t' -k1,1 -k9,9gr "$OUT/scotti_vs_swissprot.top5.tsv" \
 | awk -F '\t' '!seen[$1]++' > "$OUT/scotti_vs_swissprot.best_hit.tsv"
