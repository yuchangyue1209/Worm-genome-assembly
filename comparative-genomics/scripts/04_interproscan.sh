#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
[[ -x "${INTERPROSCAN}" ]] || { echo "ERROR: missing ${INTERPROSCAN}" >&2; exit 1; }

OUT=${FUNC}/interpro
LOG=${FUNC}/logs/interpro
mkdir -p "${OUT}" "${LOG}"

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    result=${OUT}/${sp}.tsv
    require_file "${protein}"
    [[ -s "${result}" ]] && { echo "SKIP ${sp}: ${result} exists"; continue; }
    "${INTERPROSCAN}" -i "${protein}" -b "${OUT}/${sp}" -f TSV \
        -goterms -pa -cpu "${THREADS}" > "${LOG}/${sp}.log" 2>&1
done
