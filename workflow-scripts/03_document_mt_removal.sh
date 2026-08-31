#!/usr/bin/env bash
# Document mitochondrial scaffolds removed from the assembly and reproduce the
# final nuclear FASTA. Update PRE_MT_FASTA if the immediate pre-removal file has
# a different name.
set -euo pipefail

WORK=/work/cyu/assembly/worm/hifiasm/purge_q10
PRE_MT_FASTA=${PRE_MT_FASTA:-${WORK}/worm.q10_r500k_no_contig_ec.fa}
NUCLEAR_FASTA=${WORK}/worm.q10_r500k_no_contig_ec.nuclear.final.fa
OUTDIR=${WORK}/mt-removal
REMOVED_IDS=${OUTDIR}/removed_mt_scaffolds.txt
REMOVED_FASTA=${OUTDIR}/removed_mt_scaffolds.fa
RECREATED_FASTA=${OUTDIR}/worm.q10_r500k_no_contig_ec.without_mt.fa

command -v seqkit >/dev/null || { echo "ERROR: seqkit is not in PATH" >&2; exit 1; }
for file in "${PRE_MT_FASTA}" "${NUCLEAR_FASTA}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

mkdir -p "${OUTDIR}"
tmpdir=$(mktemp -d)
trap 'rm -rf "${tmpdir}"' EXIT

grep '^>' "${PRE_MT_FASTA}" |
    sed 's/^>//; s/[[:space:]].*$//' |
    sort -u > "${tmpdir}/pre_mt.ids"

grep '^>' "${NUCLEAR_FASTA}" |
    sed 's/^>//; s/[[:space:]].*$//' |
    sort -u > "${tmpdir}/nuclear.ids"

comm -23 "${tmpdir}/pre_mt.ids" "${tmpdir}/nuclear.ids" > "${REMOVED_IDS}"
[[ -s "${REMOVED_IDS}" ]] || { echo "ERROR: no removed scaffold IDs detected" >&2; exit 1; }

seqkit grep -f "${REMOVED_IDS}" "${PRE_MT_FASTA}" > "${REMOVED_FASTA}"
seqkit grep -v -f "${REMOVED_IDS}" "${PRE_MT_FASTA}" > "${RECREATED_FASTA}"

echo "=== Removed scaffold IDs ==="
cat "${REMOVED_IDS}"
echo "=== Sequence statistics ==="
seqkit stats "${PRE_MT_FASTA}" "${REMOVED_FASTA}" "${NUCLEAR_FASTA}" "${RECREATED_FASTA}"

if cmp -s "${NUCLEAR_FASTA}" "${RECREATED_FASTA}"; then
    echo "PASS: recreated FASTA is byte-identical to the final nuclear FASTA"
else
    echo "NOTE: FASTA order/wrapping may differ; compare sequence checksums below"
    seqkit fx2tab -n -i -s "${NUCLEAR_FASTA}" | sort > "${tmpdir}/nuclear.tsv"
    seqkit fx2tab -n -i -s "${RECREATED_FASTA}" | sort > "${tmpdir}/recreated.tsv"
    diff -q "${tmpdir}/nuclear.tsv" "${tmpdir}/recreated.tsv" &&
        echo "PASS: sequence IDs and sequences are identical"
fi

echo "Removed IDs:   ${REMOVED_IDS}"
echo "Removed FASTA: ${REMOVED_FASTA}"
echo "Nuclear FASTA: ${RECREATED_FASTA}"
