#!/usr/bin/env bash
# Record the primary nuclear ncRNA scans. They are expensive and therefore run
# only when RUN_TRNASCAN=1, RUN_BARRNAP=1 or RUN_RFAM=1 is explicitly supplied.
set -euo pipefail

ASM=/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
ROOT=/work/cyu/annotation/ssol-annotation/results/ncrna_final
TRNA=${ROOT}/trna/complete
RRNA=${ROOT}/rrna
RFAM=${ROOT}/rfam_genome
DB=/work/cyu/databases/rfam
THREADS=${THREADS:-32}
mkdir -p "${TRNA}" "${RRNA}" "${RFAM}"

[[ -s "${ASM}" ]] || { echo "ERROR: missing ${ASM}" >&2; exit 1; }

if [[ "${RUN_TRNASCAN:-0}" == 1 ]]; then
    command -v tRNAscan-SE >/dev/null || { echo "ERROR: tRNAscan-SE not found" >&2; exit 1; }
    tRNAscan-SE \
        --thread "${THREADS}" \
        --eukaryotic \
        --detail \
        --output "${TRNA}/Ssolidus_tRNAscan.out" \
        --struct "${TRNA}/Ssolidus_tRNAscan.struct" \
        --bed "${TRNA}/Ssolidus_tRNAscan.bed" \
        --fasta "${TRNA}/Ssolidus_tRNA_sequences.fa" \
        "${ASM}" \
        > "${TRNA}/Ssolidus_tRNAscan.log" 2>&1
fi

if [[ "${RUN_BARRNAP:-0}" == 1 ]]; then
    command -v barrnap >/dev/null || { echo "ERROR: barrnap not found" >&2; exit 1; }
    # This installed barrnap release offers bac/arc/fun rather than an euk mode.
    barrnap --kingdom fun --threads "${THREADS}" --rrna \
        --no-trna --no-ncrna --no-mrna "${ASM}" \
        > "${RRNA}/Ssolidus_nuclear.barrnap_fun.gff3" \
        2> "${RRNA}/Ssolidus_nuclear.barrnap_fun.log"
fi

if [[ "${RUN_RFAM:-0}" == 1 ]]; then
    command -v cmscan >/dev/null || { echo "ERROR: cmscan not found" >&2; exit 1; }
    for file in "${DB}/Rfam.cm" "${DB}/Rfam.clanin"; do
        [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
    done
    /usr/bin/time -v cmscan \
        --cpu "${THREADS}" \
        --cut_ga \
        --rfam \
        --nohmmonly \
        --fmt 2 \
        --clanin "${DB}/Rfam.clanin" \
        --tblout "${RFAM}/Ssolidus_nuclear.rfam.tblout" \
        -o /dev/null \
        "${DB}/Rfam.cm" "${ASM}" \
        2> "${RFAM}/Ssolidus_nuclear.rfam.log"
    echo 0 > "${RFAM}/Ssolidus_nuclear.rfam.exitcode"
    date > "${RFAM}/Ssolidus_nuclear.rfam.finished.txt"
fi

echo "Primary scans are opt-in:"
echo "  RUN_TRNASCAN=${RUN_TRNASCAN:-0}"
echo "  RUN_BARRNAP=${RUN_BARRNAP:-0}"
echo "  RUN_RFAM=${RUN_RFAM:-0}"
echo "Use 13_ncrna_annotation_and_audit.sh to rebuild the conservative release from completed scans."
