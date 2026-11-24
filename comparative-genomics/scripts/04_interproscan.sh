#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
[[ -x "${INTERPROSCAN}" ]] || { echo "ERROR: missing ${INTERPROSCAN}" >&2; exit 1; }

OUT=${INTERPRO_OUT:-${FUNC}/interpro}
LOG=${INTERPRO_LOG:-${FUNC}/logs/interpro}
require_command python3

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    require_file "${protein}"
    python3 "${SCRIPT_DIR}/interpro_verified_run.py" \
        --protein "${protein}" --prefix "${OUT}/${sp}" \
        --executable "${INTERPROSCAN}" --threads "${THREADS}" \
        --log "${LOG}/${sp}.log"
done
