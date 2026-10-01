#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_command busco

OUT=${QC}/busco_protein
LOG=${QC}/logs/busco_protein
mkdir -p "${OUT}" "${LOG}"

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    require_file "${protein}"
    if find "${OUT}/${sp}" -maxdepth 1 -name 'short_summary*.txt' -print -quit 2>/dev/null | grep -q .; then
        echo "SKIP ${sp}: BUSCO summary exists"
        continue
    fi
    echo "START ${sp}: $(date)"
    if busco -i "${protein}" -o "${sp}" --out_path "${OUT}" \
        -m proteins -l "${BUSCO_LINEAGE}" -c "${THREADS}" \
        > "${LOG}/${sp}.log" 2>&1; then
        echo "FINISH ${sp}: $(date)"
    else
        status=$?
        echo "FAILED ${sp}: exit=${status}; see ${LOG}/${sp}.log" >&2
        tail -n 40 "${LOG}/${sp}.log" >&2
        exit "${status}"
    fi
done
