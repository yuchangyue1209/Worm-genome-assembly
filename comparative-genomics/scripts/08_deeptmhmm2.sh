#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_command dtm2

OUT=${SECRETOME}/deeptmhmm2
LOG=${SECRETOME}/logs/deeptmhmm2
BATCH_SIZE=${BATCH_SIZE:-16}
mkdir -p "${OUT}" "${LOG}"

if python -c 'import torch,sys; sys.exit(0 if torch.cuda.is_available() else 1)'; then
    DEVICE=${DEVICE:-cuda}
else
    DEVICE=${DEVICE:-cpu}
fi
echo "DeepTMHMM2 device: ${DEVICE}"

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    result=${OUT}/${sp}/predictions.json
    require_file "${protein}"
    [[ -s "${result}" ]] && { echo "SKIP ${sp}: ${result} exists"; continue; }
    mkdir -p "${OUT}/${sp}"
    dtm2 "${protein}" "${OUT}/${sp}" --device "${DEVICE}" \
        --batch-size "${BATCH_SIZE}" > "${LOG}/${sp}.log" 2>&1
    require_file "${result}"
done
