#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_file "${PROTEINS}"
[[ -x "${INTERPROSCAN}" ]] || {
    echo "ERROR: missing executable ${INTERPROSCAN}" >&2
    exit 1
}

mkdir -p "${INTERPRO_OUT}"
result=${INTERPRO_OUT}/scotti_braker.tsv
if [[ -s "${result}" ]]; then
    rows=$(wc -l < "${result}")
    echo "SKIP: validated InterProScan output exists (${rows} rows)"
    exit 0
fi

search_fasta=${INTERPRO_OUT}/scotti_braker.search.fa
python3 - "${PROTEINS}" "${search_fasta}" "${EXPECTED_PROTEINS}" <<'PY'
from pathlib import Path
import sys

source, output, expected = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
allowed = set("ACDEFGHIKLMNPQRSTVWYX")
seen = set()
records = []

with source.open() as src:
    identifier = None
    sequence = []

    def write_record() -> None:
        if identifier is None:
            return
        seq = "".join(sequence).replace(" ", "").upper().rstrip("*")
        if not seq or set(seq) - allowed:
            raise SystemExit(f"ERROR: invalid protein sequence: {identifier}")
        records.append((identifier, seq))

    for line in src:
        if line.startswith(">"):
            write_record()
            identifier = line[1:].split()[0]
            if identifier in seen:
                raise SystemExit(f"ERROR: duplicate protein ID: {identifier}")
            seen.add(identifier)
            sequence = []
        else:
            sequence.append(line.strip())
    write_record()

if len(records) != expected:
    raise SystemExit(f"ERROR: expected {expected} proteins, wrote {len(records)}")
with output.open("x") as dst:
    for identifier, seq in records:
        dst.write(f">{identifier}\n{seq}\n")
print(f"Prepared proteins: {len(records)}")
PY

"${INTERPROSCAN}" \
    -i "${search_fasta}" \
    -b "${INTERPRO_OUT}/scotti_braker" \
    -f TSV,GFF3 \
    -goterms \
    -pa \
    -cpu "${THREADS}" \
    -T "${INTERPRO_OUT}/temp" \
    2>&1 | tee "${INTERPRO_OUT}/run.log"

require_file "${result}"
echo "PASS: ${result}"
