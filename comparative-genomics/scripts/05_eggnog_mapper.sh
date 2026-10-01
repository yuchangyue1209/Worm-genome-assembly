#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_command emapper.py

for file in eggnog.db eggnog.taxa.db eggnog_proteins.dmnd; do
    require_file "${EGGNOG_DB}/${file}"
done

OUT=${FUNC}/eggnog
LOG=${FUNC}/logs/eggnog
mkdir -p "${OUT}" "${LOG}"

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    result=${OUT}/${sp}.emapper.annotations
    require_file "${protein}"
    [[ -s "${result}" ]] && { echo "SKIP ${sp}: ${result} exists"; continue; }
    emapper.py --data_dir "${EGGNOG_DB}" -i "${protein}" \
        --output "${sp}" --output_dir "${OUT}" --itype proteins \
        -m diamond --cpu "${THREADS}" --override \
        > "${LOG}/${sp}.log" 2>&1
done
