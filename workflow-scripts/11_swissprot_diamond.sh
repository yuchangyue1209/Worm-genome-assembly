#!/usr/bin/env bash
# Functional annotation of all final protein isoforms against Swiss-Prot.
# Safe to rerun only in a new/empty output directory or with OVERWRITE=1.
set -euo pipefail

PROTEINS=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_final.proteins.fa
DBFA=/work/cyu/annotation/ssol-annotation/results/functional_final/database/uniprot_sprot.fasta
ROOT=/work/cyu/annotation/ssol-annotation/results/functional_final/swissprot
DB=${ROOT}/uniprot_sprot
TOP5=${ROOT}/ssol_final_vs_swissprot.top5.tsv
BEST=${ROOT}/ssol_final_vs_swissprot.best_hit.tsv
SUMMARY=${ROOT}/swissprot_annotation_summary.tsv
THREADS=${THREADS:-32}

for cmd in diamond awk sort; do
    command -v "${cmd}" >/dev/null || { echo "ERROR: ${cmd} not found" >&2; exit 1; }
done
for file in "${PROTEINS}" "${DBFA}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done
mkdir -p "${ROOT}"

if [[ -s "${TOP5}" && "${OVERWRITE:-0}" != 1 ]]; then
    echo "ERROR: output exists: ${TOP5}; set OVERWRITE=1 to rerun" >&2
    exit 1
fi

[[ -s "${DB}.dmnd" ]] || diamond makedb --in "${DBFA}" --db "${DB}"

diamond blastp \
    --db "${DB}" \
    --query "${PROTEINS}" \
    --out "${TOP5}" \
    --threads "${THREADS}" \
    --sensitive \
    --evalue 1e-5 \
    --max-target-seqs 5 \
    --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle

# DIAMOND output is ordered by query and decreasing score, but explicit sorting
# makes the best-hit rule reproducible if output ordering changes.
sort -t$'\t' -k1,1 -k9,9gr "${TOP5}" |
    awk -F'\t' '!seen[$1]++' > "${BEST}"

total=$(grep -c '^>' "${PROTEINS}")
hits=$(awk -F'\t' 'NF&&!seen[$1]++{n++}END{print n+0}' "${BEST}")
cat > "${SUMMARY}" <<EOF
metric	value
total_final_proteins	${total}
proteins_with_swissprot_hit	${hits}
best_hit_rows	${hits}
percent_with_swissprot_hit	$(awk -v h="${hits}" -v n="${total}" 'BEGIN{printf "%.2f",100*h/n}')
EOF

column -t -s $'\t' "${SUMMARY}"
