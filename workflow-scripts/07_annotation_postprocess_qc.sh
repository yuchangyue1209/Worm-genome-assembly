#!/usr/bin/env bash
# Create the canonical BRAKER3 protein set and assess annotation completeness.
# Required tools: AGAT, BUSCO, and standard Unix utilities.
set -euo pipefail

ROOT=/work/cyu/annotation/ssol-annotation
BRAKER=${ROOT}/results/braker3_recovered
GENOME=${ROOT}/results/repeats/worm.q10_r500k_no_contig_ec.nuclear.final.fa.masked
BUSCO_OUT=${ROOT}/results/busco
THREADS=${THREADS:-32}

GFF=${BRAKER}/braker.gff3
PROTEINS=${BRAKER}/braker.aa
CANONICAL_GFF=${BRAKER}/braker.longest_isoform.gff3
CANONICAL_AA=${BRAKER}/braker.longest_isoform.aa

for cmd in agat_sp_keep_longest_isoform.pl agat_sp_extract_sequences.pl busco; do
    command -v "${cmd}" >/dev/null || {
        echo "ERROR: ${cmd} is not in PATH" >&2
        exit 1
    }
done
for file in "${GFF}" "${PROTEINS}" "${GENOME}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

if [[ ! -s "${CANONICAL_GFF}" ]]; then
    agat_sp_keep_longest_isoform.pl \
        --gff "${GFF}" \
        --output "${CANONICAL_GFF}"
fi

if [[ ! -s "${CANONICAL_AA}" ]]; then
    agat_sp_extract_sequences.pl \
        --gff "${CANONICAL_GFF}" \
        --fasta "${GENOME}" \
        --protein \
        --output "${CANONICAL_AA}"
fi

echo "genes: $(awk -F '\t' '$3=="gene"{n++} END{print n+0}' "${GFF}")"
echo "transcripts: $(awk -F '\t' '$3=="mRNA"{n++} END{print n+0}' "${GFF}")"
echo "all proteins: $(grep -c '^>' "${PROTEINS}")"
echo "canonical proteins: $(grep -c '^>' "${CANONICAL_AA}")"

awk '
    function check() {
        if (seq == "") return
        total++
        original=seq
        sub(/\*+$/, "", seq)
        if (seq ~ /\*/) internal++
        if (original !~ /\*$/) no_terminal++
        seq=""
    }
    /^>/ {check(); next}
    {gsub(/[[:space:]]/, ""); seq=seq $0}
    END {
        check()
        print "total_proteins =", total
        print "proteins_with_internal_stop =", internal+0
        print "proteins_without_terminal_stop =", no_terminal+0
    }
' "${PROTEINS}"

mkdir -p "${BUSCO_OUT}"
cd "${BUSCO_OUT}"
if [[ ! -d braker3_longest_metazoa ]]; then
    busco \
        -i "${CANONICAL_AA}" \
        -o braker3_longest_metazoa \
        -m proteins \
        -l metazoa_odb12 \
        -c "${THREADS}"
else
    echo "BUSCO output already exists: ${BUSCO_OUT}/braker3_longest_metazoa"
fi

cat "${BUSCO_OUT}"/braker3_longest_metazoa/short_summary*.txt
