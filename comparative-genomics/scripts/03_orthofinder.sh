#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_command orthofinder

INPUT=${ORTHO_ROOT}/input
mkdir -p "${INPUT}"

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    require_file "${protein}"
    ln -sfn "${protein}" "${INPUT}/${sp}.faa"
done

orthofinder -f "${INPUT}" -t "${THREADS}" -a "${THREADS}"
