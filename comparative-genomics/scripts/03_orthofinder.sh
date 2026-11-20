#!/usr/bin/env bash
# Uses a NEW frozen input set for each run; no old shared input directory is reused.
set -euo pipefail
source "$(dirname "$0")/common.sh"
require_command orthofinder
RUN_ID=${ORTHO_RUN_ID:-"background10_$(date +%Y%m%d_%H%M%S)"}
[[ "$RUN_ID" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Invalid ORTHO_RUN_ID';exit 1; }
RUNROOT="$ORTHO_ROOT/runs/$RUN_ID"
[[ ! -e "$RUNROOT" ]] || { echo "Run exists: $RUNROOT";exit 1; }
TAXA=("${SPECIES[@]}")
ADD_COTTI=${INCLUDE_SCOTTI_CANDIDATE:-0}
if [[ "$ADD_COTTI" == 1 ]];then
 for sp in "${TAXA[@]}";do [[ "$sp" != S_cotti ]] || { echo 'S_cotti already configured';exit 1; };done
 TAXA+=(S_cotti)
fi
# Validate all inputs before creating the run.
for sp in "${TAXA[@]}";do
 [[ "$sp" =~ ^[A-Za-z0-9_.-]+$ ]] || exit 1
 if [[ "$sp" == S_cotti && "$ADD_COTTI" == 1 ]];then
  STAGE=${SCOTTI_CANDIDATE_STAGE:-"${PROJECT}/02_standardized_candidates/Scotti_filtered7_v1"}
  require_file "$STAGE/STATUS.txt";protein="$STAGE/proteins_primary/S_cotti.faa"
 else protein="$STD/proteins_primary/$sp.faa";fi
 require_file "$protein"
done
mkdir -p "$RUNROOT/input"
printf 'species\tinput_source\n' > "$RUNROOT/input_manifest.tsv"
for sp in "${TAXA[@]}";do
 if [[ "$sp" == S_cotti && "$ADD_COTTI" == 1 ]];then protein="$STAGE/proteins_primary/S_cotti.faa"
 else protein="$STD/proteins_primary/$sp.faa";fi
 # Copies, rather than mutable symlinks, freeze the actual protein content for this run.
 cp "$protein" "$RUNROOT/input/$sp.faa"
 printf '%s\t%s\n' "$sp" "$protein" >> "$RUNROOT/input_manifest.tsv"
done
python3 - "$RUNROOT/input" <<'EOF'
from pathlib import Path
import sys
for p in Path(sys.argv[1]).glob('*.faa'):
 ids=set();n=0
 for line in p.open():
  if not line.startswith('>'):continue
  key=line[1:].split()[0]
  if key in ids:raise ValueError(f'Duplicate ID in {p}: {key}')
  ids.add(key);n+=1
 if not n:raise ValueError(f'No proteins: {p}')
 print(p.name,n)
EOF
sha256sum "$RUNROOT/input"/*.faa > "$RUNROOT/input.sha256"
if [[ "$ADD_COTTI" == 1 ]];then
 printf '%s\n' 'EXPLORATORY: includes cotti filtered7 candidate; not final family-count evidence.' > "$RUNROOT/STATUS.txt"
else printf '%s\n' 'BACKGROUND EXPLORATION: final inference requires curated candidate/TE review.' > "$RUNROOT/STATUS.txt";fi
orthofinder -f "$RUNROOT/input" -t "$THREADS" -a "$THREADS" 2>&1 | tee "$RUNROOT/orthofinder.log"
