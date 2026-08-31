#!/usr/bin/env bash
# Run InterProScan on all final pre-filter protein isoforms, or validate output.
set -euo pipefail

PROTEINS=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_final.proteins.fa
IPR=/work/cyu/software/interproscan-5.78-109.0
OUT=/work/cyu/annotation/ssol-annotation/results/functional_final/interproscan
PREFIX=ssol_final
THREADS=${THREADS:-32}
mkdir -p "${OUT}"

[[ -s "${PROTEINS}" ]] || { echo "ERROR: missing ${PROTEINS}" >&2; exit 1; }

if [[ "${RUN_INTERPRO:-0}" == 1 ]]; then
    [[ -x "${IPR}/interproscan.sh" ]] || { echo "ERROR: missing InterProScan" >&2; exit 1; }
    "${IPR}/interproscan.sh" \
        -i "${PROTEINS}" \
        -b "${OUT}/${PREFIX}" \
        -f TSV,GFF3 \
        -goterms \
        -pa \
        -cpu "${THREADS}"
fi

TSV=${OUT}/${PREFIX}.tsv
[[ -s "${TSV}" ]] || { echo "ERROR: missing ${TSV}" >&2; exit 1; }
any=$(awk -F'\t' 'NF&&!seen[$1]++{n++}END{print n+0}' "${TSV}")
integrated=$(awk -F'\t' '$12~/^IPR/&&!seen[$1]++{n++}END{print n+0}' "${TSV}")
printf 'any_member_database_match\t%s\nintegrated_InterPro_entry\t%s\n' "${any}" "${integrated}"
[[ "${any}" -eq 11877 && "${integrated}" -eq 10522 ]] || {
    echo "ERROR: InterPro counts differ from documented pre-filter results" >&2; exit 1;
}
echo "PASS: InterProScan annotation validated"
