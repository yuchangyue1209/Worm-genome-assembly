#!/usr/bin/env bash
# Historical batch: ten species x two lineages completed. Do not resubmit completed runs.
# Export CONFIG to the reviewed server config before using this wrapper.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
for LINEAGE in metazoa_odb10 metazoa_odb12;do
 BUSCO_LINEAGE_OVERRIDE="$LINEAGE" bash "$SCRIPT_DIR/02b_busco_proteins.sh"
done
