#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"

ACCESSION=${SPROLIFERUM_ACCESSION:-GCA_053814205.1}
SOURCE=${SPROLIFERUM_SOURCE:-"${RAW}/S_proliferum_SprJ4/ncbi_dataset/data/${ACCESSION}"}

genome=${SOURCE}/${ACCESSION}_SprJ_4.0_genomic.fna
gff=${SOURCE}/genomic.gff
cds=${SOURCE}/cds_from_genomic.fna
protein=${SOURCE}/protein.faa

for file in "${genome}" "${gff}" "${cds}" "${protein}"; do
    require_file "${file}"
done

mkdir -p "${STD}/genomes" "${STD}/gff" "${STD}/cds" \
    "${STD}/proteins_raw" "${STD}/proteins_primary"

ln -sfn "${genome}" "${STD}/genomes/S_proliferum.fasta"
ln -sfn "${gff}" "${STD}/gff/S_proliferum.gff3"
ln -sfn "${cds}" "${STD}/cds/S_proliferum.fasta"
ln -sfn "${protein}" "${STD}/proteins_raw/S_proliferum.faa"

# The public annotation has one mRNA per gene. Its protein FASTA is therefore
# already the primary protein set; untranslated models have no protein record.
ln -sfn "${protein}" "${STD}/proteins_primary/S_proliferum.faa"

protein_count=$(grep -c '^>' "${protein}")
duplicate_ids=$(awk '/^>/{id=$1; sub(/^>/,"",id); seen[id]++} END{for(id in seen) if(seen[id]>1)n++; print n+0}' "${protein}")

[[ "${protein_count}" -eq 28993 ]] || {
    echo "ERROR: expected 28,993 proteins, found ${protein_count}" >&2
    exit 1
}
[[ "${duplicate_ids}" -eq 0 ]] || {
    echo "ERROR: found ${duplicate_ids} duplicated protein IDs" >&2
    exit 1
}

echo "PASS: standardized S_proliferum (${protein_count} primary proteins)"
