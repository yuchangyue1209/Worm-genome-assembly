#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_file "${PROTEINS}"
require_file "${GTF}"
require_file "${CDS}"

protein_count=$(grep -c '^>' "${PROTEINS}")
cds_count=$(grep -c '^>' "${CDS}")
duplicate_proteins=$(awk '/^>/{id=$1; sub(/^>/,"",id); seen[id]++} END{for(id in seen) if(seen[id]>1)n++; print n+0}' "${PROTEINS}")

[[ "${protein_count}" -eq "${EXPECTED_PROTEINS}" ]] || {
    echo "ERROR: expected ${EXPECTED_PROTEINS} proteins, found ${protein_count}" >&2
    exit 1
}
[[ "${duplicate_proteins}" -eq 0 ]] || {
    echo "ERROR: found ${duplicate_proteins} duplicate protein IDs" >&2
    exit 1
}

printf 'Candidate proteins\t%s\n' "${protein_count}"
printf 'Candidate CDS sequences\t%s\n' "${cds_count}"
printf 'Duplicate protein IDs\t%s\n' "${duplicate_proteins}"
echo "PASS: S. cotti BRAKER candidate inputs validated"
