#!/usr/bin/env bash
# Integrate filtered Swiss-Prot, eggNOG, and InterPro evidence; assign reviewed
# products; and write functional attributes into the protein-coding GFF3.
set -euo pipefail

ROOT=/work/cyu/annotation/ssol-annotation
FINAL=${ROOT}/results/final_filtered_release
FUNC=${ROOT}/results/functional_final
INTEGRATION=${FINAL}/functional_integration
SCRIPTS=${ROOT}/05-functional
mkdir -p "${INTEGRATION}"

PROTEINS=${FINAL}/Ssolidus_nuclear.proteins.filtered62.fa
SWISS=${FUNC}/integration_filtered62/swissprot.filtered62.tsv
EGGNOG=${FUNC}/eggnog/ssol_final.emapper.annotations
INTERPRO=${FUNC}/integration_filtered62/interpro.filtered62.tsv
MASTER=${INTEGRATION}/Ssolidus_protein_function_master.tsv
PRODUCTS=${INTEGRATION}/Ssolidus_product_name_assignment.tsv
INPUT_GFF=${FINAL}/Ssolidus_nuclear.protein_coding.release_candidate.gff3
OUTPUT_GFF=${FINAL}/Ssolidus_nuclear.protein_coding.functional.gff3

PROGRAMS=(
    "${SCRIPTS}/05_integrate_filtered62_annotations.py"
    "${SCRIPTS}/06_assign_conservative_product_names.py"
    "${SCRIPTS}/07_write_functional_annotation_to_gff3.py"
)
INPUTS=("${PROTEINS}" "${SWISS}" "${EGGNOG}" "${INTERPRO}" "${INPUT_GFF}")
for file in "${PROGRAMS[@]}" "${INPUTS[@]}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

if [[ "${REBUILD_INTEGRATION:-0}" == 1 ]]; then
    python "${PROGRAMS[0]}" "${PROTEINS}" "${SWISS}" "${EGGNOG}" "${INTERPRO}" "${MASTER}"
    python "${PROGRAMS[1]}" "${MASTER}" "${PRODUCTS}"
    python "${PROGRAMS[2]}" "${INPUT_GFF}" "${PRODUCTS}" "${MASTER}" "${OUTPUT_GFF}"
fi

for file in "${MASTER}" "${PRODUCTS}" "${OUTPUT_GFF}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing completed output ${file}" >&2; exit 1; }
done

master_n=$(awk 'END{print NR-1}' "${MASTER}")
product_n=$(awk 'END{print NR-1}' "${PRODUCTS}")
assigned=$(awk -F'\t' 'NR>1&&$4!="hypothetical protein"{n++}END{print n+0}' "${PRODUCTS}")
hypothetical=$(awk -F'\t' 'NR>1&&$4=="hypothetical protein"{n++}END{print n+0}' "${PRODUCTS}")
manual=$(awk -F'\t' 'NR>1&&$14=="yes"{n++}END{print n+0}' "${PRODUCTS}")
transcript_products=$(awk -F'\t' '$3=="transcript"&&$9~/(^|;)product=/{n++}END{print n+0}' "${OUTPUT_GFF}")
cds_products=$(awk -F'\t' '$3=="CDS"&&$9~/(^|;)product=/{n++}END{print n+0}' "${OUTPUT_GFF}")

[[ "${master_n}" -eq 12849 && "${product_n}" -eq 12849 ]] || { echo "ERROR: integration row count" >&2; exit 1; }
[[ "${assigned}" -eq 10535 && "${hypothetical}" -eq 2314 && "${manual}" -eq 0 ]] || { echo "ERROR: product QC" >&2; exit 1; }
[[ "${transcript_products}" -eq 12849 && "${cds_products}" -eq 104603 ]] || { echo "ERROR: GFF3 product coverage" >&2; exit 1; }

printf 'integrated_proteins\t%s\nassigned_products\t%s\nhypothetical_proteins\t%s\nmanual_review_remaining\t%s\n' \
    "${master_n}" "${assigned}" "${hypothetical}" "${manual}"
echo "PASS: functional integration validated"
