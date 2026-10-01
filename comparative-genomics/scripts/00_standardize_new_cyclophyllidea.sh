#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"

mkdir -p "${STD}"/{genomes,gff,cds,proteins_raw,proteins_primary}

eg=${RAW}/E_granulosus/ncbi_dataset/data/GCA_021556725.1
ts=${RAW}/T_solium/ncbi_dataset/data/GCA_001870725.1
tm=${RAW}/T_multiceps/WBPS19

for file in \
    "${eg}/GCA_021556725.1_ASM2155672v1_genomic.fna" \
    "${eg}/genomic.gff" "${eg}/cds_from_genomic.fna" "${eg}/protein.faa" \
    "${ts}/GCA_001870725.1_MEX_genome_complete.1-6-13_genomic.fna" \
    "${ts}/genomic.gff" "${ts}/cds_from_genomic.fna" "${ts}/protein.faa" \
    "${tm}/taenia_multiceps.PRJNA307624.WBPS19.genomic.fa.gz" \
    "${tm}/taenia_multiceps.PRJNA307624.WBPS19.annotations.gff3.gz" \
    "${tm}/taenia_multiceps.PRJNA307624.WBPS19.CDS_transcripts.fa.gz" \
    "${tm}/taenia_multiceps.PRJNA307624.WBPS19.protein.fa.gz"
do
    require_file "${file}"
done

ln -sfn "${eg}/GCA_021556725.1_ASM2155672v1_genomic.fna" "${STD}/genomes/E_granulosus.fasta"
ln -sfn "${eg}/genomic.gff" "${STD}/gff/E_granulosus.gff3"
ln -sfn "${eg}/cds_from_genomic.fna" "${STD}/cds/E_granulosus.fasta"
ln -sfn "${eg}/protein.faa" "${STD}/proteins_raw/E_granulosus.faa"

ln -sfn "${ts}/GCA_001870725.1_MEX_genome_complete.1-6-13_genomic.fna" "${STD}/genomes/T_solium.fasta"
ln -sfn "${ts}/genomic.gff" "${STD}/gff/T_solium.gff3"
ln -sfn "${ts}/cds_from_genomic.fna" "${STD}/cds/T_solium.fasta"
ln -sfn "${ts}/protein.faa" "${STD}/proteins_raw/T_solium.faa"

gzip -dc "${tm}/taenia_multiceps.PRJNA307624.WBPS19.genomic.fa.gz" > "${STD}/genomes/T_multiceps.fasta.tmp"
gzip -dc "${tm}/taenia_multiceps.PRJNA307624.WBPS19.annotations.gff3.gz" > "${STD}/gff/T_multiceps.gff3.tmp"
gzip -dc "${tm}/taenia_multiceps.PRJNA307624.WBPS19.CDS_transcripts.fa.gz" > "${STD}/cds/T_multiceps.fasta.tmp"
gzip -dc "${tm}/taenia_multiceps.PRJNA307624.WBPS19.protein.fa.gz" > "${STD}/proteins_raw/T_multiceps.faa.tmp"
mv "${STD}/genomes/T_multiceps.fasta.tmp" "${STD}/genomes/T_multiceps.fasta"
mv "${STD}/gff/T_multiceps.gff3.tmp" "${STD}/gff/T_multiceps.gff3"
mv "${STD}/cds/T_multiceps.fasta.tmp" "${STD}/cds/T_multiceps.fasta"
mv "${STD}/proteins_raw/T_multiceps.faa.tmp" "${STD}/proteins_raw/T_multiceps.faa"

# These releases contain one protein per protein-coding transcript/gene.
for sp in E_granulosus T_solium T_multiceps; do
    ln -sfn "${STD}/proteins_raw/${sp}.faa" "${STD}/proteins_primary/${sp}.faa"
done

echo "PASS: standardized E_granulosus, T_solium, and T_multiceps"
