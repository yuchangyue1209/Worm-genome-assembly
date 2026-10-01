#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_command busco

OUT=${QC}/busco_genome
mkdir -p "${OUT}"

for sp in "${SPECIES[@]}"; do
    genome=${STD}/genomes/${sp}.fasta
    require_file "${genome}"
    if find "${OUT}/${sp}" -name 'short_summary*.txt' -print -quit 2>/dev/null | grep -q .; then
        echo "SKIP ${sp}: BUSCO summary exists"
        continue
    fi
    busco -i "${genome}" -o "${sp}" --out_path "${OUT}" \
        -m genome -l "${BUSCO_LINEAGE}" -c "${THREADS}"
done
