#!/usr/bin/env bash
# De novo repeat discovery and soft masking of the final nuclear assembly.
set -euo pipefail

ROOT=/work/cyu/annotation/ssol-annotation
GENOME=/work/cyu/assembly/worm/hifiasm/purge_q10/worm.q10_r500k_no_contig_ec.nuclear.final.fa
LABEL=ssol
THREADS=32
OUT=${ROOT}/results/repeats
TMP=${ROOT}/temp/repeatmodeler
mkdir -p "${OUT}" "${TMP}" "${ROOT}/logs"

for cmd in BuildDatabase RepeatModeler RepeatMasker; do
    command -v "${cmd}" >/dev/null || { echo "ERROR: ${cmd} is not in PATH" >&2; exit 1; }
done
[[ -s "${GENOME}" ]] || { echo "ERROR: missing ${GENOME}" >&2; exit 1; }

BuildDatabase -name "${OUT}/${LABEL}" "${GENOME}"
cd "${TMP}"
RepeatModeler -database "${OUT}/${LABEL}" -threads "${THREADS}" -LTRStruct \
    2>&1 | tee "${ROOT}/logs/repeatmodeler.log"

# RepeatModeler normally writes this beside the database. If classification was
# completed manually, copy consensi.fa.classified here before continuing.
LIB=${OUT}/${LABEL}-families.fa
if [[ ! -s "${LIB}" ]]; then
    classified=$(find "${TMP}" -type f -name 'consensi.fa.classified' -print | head -n 1)
    [[ -n "${classified}" ]] || { echo "ERROR: classified repeat library not found" >&2; exit 1; }
    cp "${classified}" "${LIB}"
fi

RepeatMasker -pa "${THREADS}" -lib "${LIB}" -xsmall -gff \
    -dir "${OUT}" "${GENOME}" \
    2>&1 | tee "${ROOT}/logs/repeatmasker.log"

base=$(basename "${GENOME}")
ls -lh \
    "${OUT}/${base}.masked" \
    "${OUT}/${base}.out" \
    "${OUT}/${base}.out.gff" \
    "${OUT}/${base}.tbl"
