#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_file "${PROTEINS}"

printf 'Dataset\tStatus\tUnique_proteins_or_rows\n'
printf 'BRAKER_proteins\tCOMPLETE\t%s\n' "$(grep -c '^>' "${PROTEINS}")"

eggnog=${EGGNOG_OUT}/scotti_braker.emapper.annotations
if [[ -s "${eggnog}" ]]; then
    count=$(awk '$0 !~ /^#/ && NF && !seen[$1]++ {n++} END{print n+0}' "${eggnog}")
    printf 'eggNOG\tCOMPLETE\t%s\n' "${count}"
else
    printf 'eggNOG\tMISSING\t0\n'
fi

interpro=${INTERPRO_OUT}/scotti_braker.tsv
if [[ -s "${interpro}" ]]; then
    count=$(awk -F '\t' 'NF && !seen[$1]++ {n++} END{print n+0}' "${interpro}")
    printf 'InterProScan\tCOMPLETE\t%s\n' "${count}"
else
    printf 'InterProScan\tMISSING_OR_RUNNING\t0\n'
fi
