#!/usr/bin/env bash
# Audit the final proteome and summarize the completed Swiss-Prot annotation.
set -euo pipefail

PROTEINS=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_final.proteins.fa
SWISS=/work/cyu/annotation/ssol-annotation/results/functional_final/swissprot/ssol_final_vs_swissprot.best_hit.tsv
OUT=/work/cyu/annotation/ssol-annotation/results/release_qc/02_protein_audit
mkdir -p "${OUT}"

for file in "${PROTEINS}" "${SWISS}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

awk -v out="${OUT}" -v summary="${OUT}/protein_audit_summary.tsv" '
function finish(){
    if(id=="") return
    n++; len=length(seq); total+=len
    print id,len > lengths
    if(len<50) print id > short_ids
    if(len>5000) print id > long_ids
    x=seq; sub(/\*+$/, "", x); if(x ~ /\*/) internal++
}
BEGIN {
    lengths=out "/protein_lengths.tsv"
    short_ids=out "/proteins_lt50aa.txt"
    long_ids=out "/proteins_gt5000aa.txt"
}
/^>/ {finish(); id=substr($1,2); seq=""; next}
{gsub(/[[:space:]]/,""); seq=seq $0}
END {
    finish()
    print "metric\tvalue" > summary
    print "total_proteins\t" n >> summary
    print "total_amino_acids\t" total >> summary
    print "mean_length_aa\t" total/n >> summary
    print "proteins_with_internal_stop\t" internal+0 >> summary
}
' "${PROTEINS}"

printf 'proteins_lt50aa\t%s\n' "$(wc -l < "${OUT}/proteins_lt50aa.txt")" >> "${OUT}/protein_audit_summary.tsv"
printf 'proteins_gt5000aa\t%s\n' "$(wc -l < "${OUT}/proteins_gt5000aa.txt")" >> "${OUT}/protein_audit_summary.tsv"

awk -F'\t' 'NR==FNR{wanted[$1]=1; next} ($1 in wanted){print}' \
    "${OUT}/proteins_lt50aa.txt" "${SWISS}" > "${OUT}/proteins_lt50aa.with_swissprot.tsv"
awk -F'\t' 'NR==FNR{wanted[$1]=1; next} ($1 in wanted){print}' \
    "${OUT}/proteins_gt5000aa.txt" "${SWISS}" > "${OUT}/proteins_gt5000aa.with_swissprot.tsv"

total=$(grep -c '^>' "${PROTEINS}")
hits=$(awk -F'\t' 'NF && !seen[$1]++{n++} END{print n+0}' "${SWISS}")
printf 'proteins_with_swissprot_hit\t%s\npercent_with_swissprot_hit\t%.2f\n' \
    "${hits}" "$(awk -v h="${hits}" -v n="${total}" 'BEGIN{print 100*h/n}')" \
    >> "${OUT}/protein_audit_summary.tsv"

column -t -s $'\t' "${OUT}/protein_audit_summary.tsv"
