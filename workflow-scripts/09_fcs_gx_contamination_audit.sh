#!/usr/bin/env bash
# Reproduce or document the NCBI FCS-GX screen of the final nuclear assembly.
# RUN_FCS=1 performs the expensive database scan; the default only validates
# existing results and the cleaned release.
set -euo pipefail

ASM=/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
GFF=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3
GTF=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.checked.gtf
DB=/work/cyu/databases/ncbi_fcs_gx/gxdb/all
OUT=/work/cyu/annotation/ssol-annotation/results/release_qc/01_contamination/fcs_gx
REPORT=${OUT}/Ssolidus_nuclear.70667.fcs_gx_report.txt

mkdir -p "${OUT}"

if [[ "${RUN_FCS:-0}" == 1 ]]; then
    command -v run_gx.py >/dev/null || { echo "ERROR: run_gx.py not found" >&2; exit 1; }
    [[ -s "${DB}.gxi" ]] || { echo "ERROR: missing FCS-GX database ${DB}.gxi" >&2; exit 1; }
    run_gx.py --fasta="${ASM}" --tax-id=70667 --gx-db="${DB}" --out-dir="${OUT}"
fi

[[ -s "${REPORT}" ]] || { echo "ERROR: missing existing FCS report ${REPORT}" >&2; exit 1; }
echo "=== FCS-GX actions ==="
grep -v '^#' "${REPORT}"

echo "=== Verify excluded scaffold has no annotations ==="
for file in "${GFF}" "${GTF}"; do
    n=$(awk -F'\t' '$1=="scaffold_61"{n++} END{print n+0}' "${file}")
    printf '%s\t%s records\n' "${file}" "${n}"
    [[ "${n}" -eq 0 ]] || { echo "ERROR: scaffold_61 still has annotations" >&2; exit 1; }
done

if grep -q '^>scaffold_61\([[:space:]]\|$\)' "${ASM}"; then
    echo "ERROR: scaffold_61 is still present in the release FASTA" >&2
    exit 1
fi

sha256sum "${ASM}"
echo "PASS: FCS contaminant scaffold_61 is absent and had no annotations"
