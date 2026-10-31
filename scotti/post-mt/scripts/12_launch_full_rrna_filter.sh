#!/usr/bin/env bash
# Environment: rrna_filter. Only launch if no full run is already active.
# Use /work because /scratch/cyu returned EIO on 2026-10-28.
source "$(dirname "$0")/00_config.sh"
# In the historical attempt these two directories were confirmed EMPTY.
# Run rmdir manually ONLY if both are still empty; never replace with rm -rf:
# rmdir "$ANNO/results/rnaseq_qc/05_sortmerna_full"
# rmdir "$ANNO/data/rnaseq/rrna_screened_v1"
python3 "$(dirname "$0")/11_run_rrna_chunks.py" \
  --scratch /work/cyu/scotti_sortmerna_tmp_v1 --chunk-pairs 2000000 \
  2>&1 | tee "$ANNO/sortmerna_full.work_tmp.console.log"
# Resume requires existing configuration.json and the SAME inputs/chunk size:
# python3 11_run_rrna_chunks.py --scratch /work/cyu/scotti_sortmerna_tmp_v1 --chunk-pairs 2000000 --resume
