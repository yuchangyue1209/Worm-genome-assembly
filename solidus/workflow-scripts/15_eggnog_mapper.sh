#!/usr/bin/env bash
# Run eggNOG-mapper on the final pre-filter protein set, or validate the result.
set -euo pipefail

PROTEINS=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_final.proteins.fa
DB=/work/cyu/databases/eggnog_mapper
OUT=/work/cyu/annotation/ssol-annotation/results/functional_final/eggnog
PREFIX=ssol_final
THREADS=${THREADS:-32}
mkdir -p "${OUT}"

for file in "${PROTEINS}" "${DB}/eggnog.db" "${DB}/eggnog.taxa.db" "${DB}/eggnog_proteins.dmnd"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

if [[ "${RUN_EGGNOG:-0}" == 1 ]]; then
    command -v emapper.py >/dev/null || { echo "ERROR: emapper.py not found" >&2; exit 1; }
    emapper.py \
        --data_dir "${DB}" \
        --input "${PROTEINS}" \
        --output "${PREFIX}" \
        --output_dir "${OUT}" \
        --itype proteins \
        --method diamond \
        --cpu "${THREADS}" \
        --override
fi

ANNOT=${OUT}/${PREFIX}.emapper.annotations
[[ -s "${ANNOT}" ]] || { echo "ERROR: missing ${ANNOT}" >&2; exit 1; }
rows=$(awk '$0!~/^#/{n++}END{print n+0}' "${ANNOT}")
[[ "${rows}" -eq 8177 ]] || { echo "ERROR: expected 8177 eggNOG rows, found ${rows}" >&2; exit 1; }
printf 'eggNOG_annotation_rows\t%s\n' "${rows}"
echo "PASS: eggNOG annotation validated"
