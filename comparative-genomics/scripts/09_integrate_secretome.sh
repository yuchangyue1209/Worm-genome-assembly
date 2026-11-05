#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"

TM_START_CUTOFF=${TM_START_CUTOFF:-40}
OUT=${SECRETOME}/integrated_final

python "${SCRIPT_DIR}/09_integrate_secretome.py" \
    --protein-dir "${STD}/proteins_primary" \
    --signalp-dir "${SECRETOME}/signalp6/results" \
    --deeptm-dir "${SECRETOME}/deeptmhmm2" \
    --output-dir "${OUT}" \
    --tm-start-cutoff "${TM_START_CUTOFF}" \
    --species "${SPECIES[@]}"

column -t -s $'\t' "${OUT}/secretome_summary.tsv"
