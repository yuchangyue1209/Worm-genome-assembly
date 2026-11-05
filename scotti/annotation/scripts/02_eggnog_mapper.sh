#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_command emapper.py
require_file "${PROTEINS}"

for file in eggnog.db eggnog.taxa.db eggnog_proteins.dmnd; do
    require_file "${EGGNOG_DB}/${file}"
done

mkdir -p "${EGGNOG_OUT}"
result=${EGGNOG_OUT}/scotti_braker.emapper.annotations
if [[ -s "${result}" ]]; then
    rows=$(awk '$0 !~ /^#/ && NF {n++} END{print n+0}' "${result}")
    echo "SKIP: validated eggNOG output exists (${rows} annotation rows)"
    exit 0
fi

emapper.py \
    --data_dir "${EGGNOG_DB}" \
    -i "${PROTEINS}" \
    --output scotti_braker \
    --output_dir "${EGGNOG_OUT}" \
    --itype proteins \
    -m diamond \
    --cpu "${THREADS}" \
    --override \
    > "${EGGNOG_OUT}/run.log" 2>&1

require_file "${result}"
echo "PASS: ${result}"
