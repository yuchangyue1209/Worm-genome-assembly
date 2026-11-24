#!/usr/bin/env bash
# Freeze a NEW run; never relink the input of a completed analysis.
# Uses the configured species set, with optional explicit release/candidate input.
set -euo pipefail
source "$(dirname "$0")/common.sh"
require_command orthofinder
require_command sha256sum
RUN_ID=${ORTHO_RUN_ID:-"configured_baseline_$(date +%Y%m%d_%H%M%S)"}
[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid ORTHO_RUN_ID'; exit 1; }
RUNROOT="$ORTHO_ROOT/runs/$RUN_ID"
[[ ! -e "$RUNROOT" ]] || { echo "Run exists: $RUNROOT"; exit 1; }

ADD_COTTI=${INCLUDE_SCOTTI_CANDIDATE:-0}
ADD_RELEASE=${INCLUDE_SCOTTI_RELEASE:-0}
[[ "$ADD_COTTI" =~ ^[01]$ && "$ADD_RELEASE" =~ ^[01]$ ]] || {
    echo 'Cotti inclusion flags must be 0 or 1'; exit 1;
}
[[ "$ADD_COTTI$ADD_RELEASE" != 11 ]] || {
    echo 'Choose release mode or historical candidate mode, not both'; exit 1;
}
TAXA=("${SPECIES[@]}")
if [[ "$ADD_COTTI" == 1 || "$ADD_RELEASE" == 1 ]]; then
    found_cotti=0
    for sp in "${TAXA[@]}"; do
        [[ "$sp" != S_cotti ]] || found_cotti=1
    done
    if [[ "$found_cotti" == 0 ]]; then
        TAXA+=(S_cotti)
    fi
fi

if [[ "$ADD_RELEASE" == 1 ]]; then
    RELEASE=${SCOTTI_RELEASE:-/work/shared/Scotti_release_v1.0_preNCBI_20261118}
    require_file "$RELEASE/summary.json"
    require_file "$RELEASE/SHA256SUMS"
    # Check the selected published baseline and the precise primary-file checksum.
    python3 - "$RELEASE" <<'PY'
from pathlib import Path
import hashlib,json,sys
root=Path(sys.argv[1])
s=json.loads((root/'summary.json').read_text())
if s.get('status')!='scientific_release_v1.0_preNCBI':
    raise ValueError('Expected the published scientific v1.0 release')
if (s.get('genes'),s.get('all_isoform_proteins'),s.get('primary_proteins'))!=(9234,11464,9234):
    raise ValueError('Unexpected selected release protein counts')
name='proteins/Scotti_proteins.longest_isoform.v1.0.fa'
entries=[line.rstrip('\n').split('  ',1)[0] for line in (root/'SHA256SUMS').open()
         if '  ' in line and line.rstrip('\n').split('  ',1)[1]==name]
if len(entries)!=1:
    raise ValueError('Missing/ambiguous release primary-protein checksum')
h=hashlib.sha256()
with (root/name).open('rb') as f:
    for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
if h.hexdigest()!=entries[0]:
    raise ValueError('Release primary FASTA checksum mismatch')
PY
fi

INPUTS=()
CHECK_ARGS=()
for sp in "${TAXA[@]}"; do
    [[ "$sp" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Invalid species: $sp"; exit 1; }
    if [[ "$sp" == S_cotti && "$ADD_RELEASE" == 1 ]]; then
        protein="$RELEASE/proteins/Scotti_proteins.longest_isoform.v1.0.fa"
    elif [[ "$sp" == S_cotti && "$ADD_COTTI" == 1 ]]; then
        STAGE=${SCOTTI_CANDIDATE_STAGE:-"${PROJECT}/02_standardized_candidates/Scotti_filtered7_v1"}
        require_file "$STAGE/STATUS.txt"
        protein="$STAGE/proteins_primary/S_cotti.faa"
    else
        protein="$STD/proteins_primary/$sp.faa"
    fi
    require_file "$protein"
    INPUTS+=("$protein")
    CHECK_ARGS+=("$sp" "$protein")
done
# Reject duplicate taxa, duplicate IDs and empty proteins before allocating a run.
python3 - "${CHECK_ARGS[@]}" <<'PY'
from pathlib import Path
import hashlib,sys
args=sys.argv[1:]
taxa=set()
for sp,name in zip(args[::2],args[1::2]):
    if sp in taxa:raise ValueError('Duplicate taxon: '+sp)
    taxa.add(sp);ids=set();current=None;length=0
    with Path(name).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if current is not None and not length:raise ValueError('Empty protein: '+current)
                tokens=line[1:].split()
                if not tokens:raise ValueError('Empty FASTA header: '+name)
                current=tokens[0];length=0
                if current in ids:raise ValueError('Duplicate protein ID: '+current)
                ids.add(current)
            elif line.strip():
                if current is None:raise ValueError('Sequence before header: '+name)
                length+=len(''.join(line.split()))
    if not ids or not length:raise ValueError('Empty FASTA/last protein: '+name)
    if sp=='S_cotti':
        if len(ids)!=9234:raise ValueError('Baseline requires 9234 cotti representatives')
        digest=hashlib.sha256(Path(name).read_bytes()).hexdigest()
        if digest!='07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371':
            raise ValueError('S_cotti input differs from the selected 9234-protein release')
    print(sp,len(ids))
PY

# Reserve this run atomically; refuse existing or concurrently created runs.
mkdir -p "$ORTHO_ROOT/runs"
mkdir "$RUNROOT"
mkdir "$RUNROOT/input"
printf 'species\tinput_source\n' > "$RUNROOT/input_manifest.tsv"
for i in "${!TAXA[@]}"; do
    sp=${TAXA[$i]}
    protein=${INPUTS[$i]}
    cp "$protein" "$RUNROOT/input/$sp.faa"
    source_sha=$(sha256sum "$protein" | awk '{print $1}')
    copied_sha=$(sha256sum "$RUNROOT/input/$sp.faa" | awk '{print $1}')
    [[ "$source_sha" == "$copied_sha" ]] || { echo "Copy checksum mismatch: $sp"; exit 1; }
    printf '%s\t%s\n' "$sp" "$protein" >> "$RUNROOT/input_manifest.tsv"
done
sha256sum "$RUNROOT/input"/*.faa > "$RUNROOT/input.sha256"
if [[ "$ADD_RELEASE" == 1 ]]; then
    printf '%s\n' \
        'BASELINE: includes the published cotti v1.0 release, 9234 representatives; all 32 TE candidates retained.' \
        'Comparative interpretation still requires cross-species annotation/family review. No uniform TE filtering is implied.' \
        > "$RUNROOT/STATUS.txt"
elif [[ "$ADD_COTTI" == 1 ]]; then
    printf '%s\n' 'HISTORICAL CANDIDATE INPUT: includes cotti filtered7; not a uniform TE-filtered analysis.' > "$RUNROOT/STATUS.txt"
elif [[ " ${TAXA[*]} " == *" S_cotti "* ]]; then
    printf '%s\n' \
        'BASELINE: configured S_cotti input matches the selected 9234-protein release checksum; all 32 TE candidates retained.' \
        'No uniform TE filtering is implied.' > "$RUNROOT/STATUS.txt"
else
    printf '%s\n' 'BACKGROUND EXPLORATION: S_cotti is not in this configured species set.' > "$RUNROOT/STATUS.txt"
fi
# Preserve the existing search/analysis thread parameters. No restart flags are inferred.
orthofinder -f "$RUNROOT/input" -t "$THREADS" -a "$THREADS" 2>&1 | tee "$RUNROOT/orthofinder.log"
