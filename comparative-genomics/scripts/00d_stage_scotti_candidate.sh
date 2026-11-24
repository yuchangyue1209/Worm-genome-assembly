#!/usr/bin/env bash
# Link the selected 9,234-protein S. cotti pre-NCBI scientific release into the
# standardized comparative-genomics layout. NCBI accessioning is a separate
# release step.
set -euo pipefail

ROOT=${ROOT:-/work/cyu/Schistocephalus_comparative_genomics}
STD=${STD:-"${ROOT}/02_standardized_data"}
COTTI_RELEASE_ROOT=${COTTI_RELEASE_ROOT:-/work/shared/Scotti_release_v1.0_preNCBI_20261118}
COTTI_GENOME=${COTTI_GENOME:-${COTTI_RELEASE_ROOT}/assembly/Scotti_nuclear.fa}
COTTI_GFF=${COTTI_GFF:-${COTTI_RELEASE_ROOT}/annotation/Scotti_nuclear.protein_coding.v1.0.gff3}
COTTI_PROTEINS=${COTTI_PROTEINS:-${COTTI_RELEASE_ROOT}/proteins/Scotti_proteins.longest_isoform.v1.0.fa}
COTTI_ALL_PROTEINS=${COTTI_ALL_PROTEINS:-${COTTI_RELEASE_ROOT}/proteins/Scotti_proteins.all_isoforms.v1.0.fa}
COTTI_CDS=${COTTI_CDS:-}
EXPECTED_PROTEINS=${EXPECTED_PROTEINS:-9234}
EXPECTED_ALL_PROTEINS=${EXPECTED_ALL_PROTEINS:-11464}
EXPECTED_SHA256=${EXPECTED_SHA256:-07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371}

for file in "$COTTI_GENOME" "$COTTI_GFF" "$COTTI_PROTEINS"; do
    [[ -s "$file" ]] || { echo "ERROR: missing or empty $file" >&2; exit 1; }
done

count=$(awk '/^>/{n++} END{print n+0}' "$COTTI_PROTEINS")
duplicates=$(awk '/^>/{id=$1;sub(/^>/,"",id);seen[id]++} END{for(id in seen)if(seen[id]>1)n++;print n+0}' "$COTTI_PROTEINS")
[[ "$count" -eq "$EXPECTED_PROTEINS" ]] || {
    echo "ERROR: expected $EXPECTED_PROTEINS S. cotti proteins; found $count" >&2
    exit 1
}
[[ "$duplicates" -eq 0 ]] || {
    echo "ERROR: S. cotti protein FASTA contains $duplicates duplicate IDs" >&2
    exit 1
}
observed_sha256=$(sha256sum "$COTTI_PROTEINS" | awk '{print $1}')
[[ "$observed_sha256" == "$EXPECTED_SHA256" ]] || {
    echo "ERROR: unexpected S. cotti release SHA256: $observed_sha256" >&2
    exit 1
}

# Preflight all inputs and destinations before creating links; never overwrite.
python3 - "$STD" "$COTTI_GENOME" "$COTTI_GFF" "$COTTI_PROTEINS" \
    "$COTTI_ALL_PROTEINS" "$COTTI_CDS" "$EXPECTED_ALL_PROTEINS" <<'PY'
import hashlib
from pathlib import Path
import sys

std, genome, gff, primary, raw, cds, expected = sys.argv[1:]
std = Path(std)

def fasta(path):
    records = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            fields = line[1:].split()
            if not fields or fields[0] in records:
                raise ValueError(f'Missing/duplicate protein ID: {line} in {path}')
            name = fields[0]
            records[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError(f'Sequence before header: {path}')
            records[name] += line.strip()
    if not records or any(not value for value in records.values()):
        raise ValueError(f'Empty protein sequence/set: {path}')
    return records

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

all_proteins = fasta(raw)
representatives = fasta(primary)
if len(all_proteins) != int(expected):
    raise ValueError(f'Expected {expected} all-isoform proteins; found {len(all_proteins)}')
if any(all_proteins.get(key) != sequence for key, sequence in representatives.items()):
    raise ValueError('Representative IDs/sequences must match the all-isoform release')
pairs = [(genome, 'genomes/S_cotti.fasta'), (gff, 'gff/S_cotti.gff3'),
         (primary, 'proteins_primary/S_cotti.faa'), (raw, 'proteins_raw/S_cotti.faa')]
if cds:
    pairs.append((cds, 'cds/S_cotti.fasta'))
for source, relative in pairs:
    source, target = Path(source), std/relative
    if not source.is_file() or not source.stat().st_size:
        raise ValueError(f'Missing or empty input: {source}')
    if target.exists() or target.is_symlink():
        if not target.is_file() or digest(target) != digest(source):
            raise ValueError(f'Conflicting destination; not overwritten: {target}')
for source, relative in pairs:
    target = std/relative
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(Path(source).resolve())
    print('VERIFIED:', target)
if not cds:
    print('NOTE: COTTI_CDS unset; existing CDS data left untouched')
PY

for target in \
    "$STD/genomes/S_cotti.fasta" \
    "$STD/gff/S_cotti.gff3" \
    "$STD/proteins_primary/S_cotti.faa"; do
    [[ -s "$target" ]] || { echo "ERROR: staged link is invalid: $target" >&2; exit 1; }
done

printf 'Species\tProteins\tRelease_version\tSHA256\tStatus\n'
printf 'S_cotti\t%s\tScotti_release_v1.0_preNCBI_20261118\t%s\tselected_scientific_release_baseline\n' \
    "$count" "$observed_sha256"
