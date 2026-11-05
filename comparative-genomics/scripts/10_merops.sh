#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
MEROPS_DB=${MEROPS_DB:-/work/cyu/databases/merops_12.4/merops_pepunit}
MEROPS_FASTA=${MEROPS_FASTA:-/work/cyu/databases/merops_12.4/pepunit.lib}
require_command blastp
require_file "${MEROPS_FASTA}"
[[ -s "${MEROPS_DB}.pin" ]] || {
    echo "ERROR: missing BLAST database ${MEROPS_DB}.pin" >&2
    exit 1
}

OUT=${SECRETOME}/merops
SECRETOME_INPUT=${SECRETOME}/integrated_final
RAW_OUT=${OUT}/raw_classical_secretome
LOG=${OUT}/logs_classical_secretome
EVALUE=${MEROPS_EVALUE:-1e-5}
SUBJECT_COVERAGE=${MEROPS_SUBJECT_COVERAGE:-50}
MAX_TARGET_SEQS=${MEROPS_MAX_TARGET_SEQS:-50}
mkdir -p "${RAW_OUT}" "${LOG}"

for sp in "${SPECIES[@]}"; do
    protein=${SECRETOME_INPUT}/${sp}.classical_secretome.faa
    result=${RAW_OUT}/${sp}.merops.blastp.tsv
    require_file "${protein}"
    if [[ -s "${result}" ]]; then
        echo "SKIP ${sp}: ${result} exists"
        continue
    fi
    echo "START ${sp}: $(date)"
    blastp \
        -query "${protein}" \
        -db "${MEROPS_DB}" \
        -evalue "${EVALUE}" \
        -max_target_seqs "${MAX_TARGET_SEQS}" \
        -seg yes \
        -num_threads "${THREADS}" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
        -out "${result}" \
        > "${LOG}/${sp}.log" 2>&1
    require_file "${result}"
    echo "FINISH ${sp}: $(date), hits=$(wc -l < "${result}")"
done

python "${SCRIPT_DIR}/10_parse_merops.py" \
    --raw-dir "${RAW_OUT}" \
    --merops-fasta "${MEROPS_FASTA}" \
    --output "${OUT}/merops_confirmed_secretome_annotations_all_species.tsv" \
    --summary "${OUT}/merops_confirmed_secretome_summary_all_species.tsv" \
    --evalue "${EVALUE}" \
    --subject-coverage "${SUBJECT_COVERAGE}" \
    --species "${SPECIES[@]}"

column -t -s $'\t' "${OUT}/merops_confirmed_secretome_summary_all_species.tsv"
