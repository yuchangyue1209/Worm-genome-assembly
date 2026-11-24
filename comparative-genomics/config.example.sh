#!/usr/bin/env bash

# Copy to config.sh and edit for the server. Do not commit private paths or data.
PROJECT=${PROJECT:-/work/cyu/Schistocephalus_comparative_genomics}
RAW=${RAW:-"${PROJECT}/01_raw_data"}
STD=${STD:-"${PROJECT}/02_standardized_data"}
QC=${QC:-"${PROJECT}/03_qc"}
ORTHO_ROOT=${ORTHO_ROOT:-"${PROJECT}/04_orthofinder"}
FUNC=${FUNC:-"${PROJECT}/05_functional_annotation"}
FAMILY_QC=${FAMILY_QC:-"${PROJECT}/05_gene_family_qc"}
SECRETOME=${SECRETOME:-"${PROJECT}/06_secretome"}

INTERPROSCAN=${INTERPROSCAN:-/work/cyu/software/interproscan-5.78-109.0/interproscan.sh}
EGGNOG_DB=${EGGNOG_DB:-/work/cyu/databases/eggnog_mapper}
MEROPS_DB=${MEROPS_DB:-/work/cyu/databases/merops_12.4/merops_pepunit}
MEROPS_FASTA=${MEROPS_FASTA:-/work/cyu/databases/merops_12.4/pepunit.lib}
BUSCO_LINEAGE=${BUSCO_LINEAGE:-metazoa_odb12}
THREADS=${THREADS:-32}

# odb12 is the primary protein-completeness comparison; odb10 is retained as
# a historical cross-version sensitivity analysis. Never combine the two
# lineage datasets in one numerical table.
BUSCO_LINEAGES=(metazoa_odb10 metazoa_odb12)

# Current eleven-species analysis set. S_cotti uses the selected 9,234-protein
# pre-NCBI scientific release. NCBI accessioning remains a separate step.
SPECIES=(
  S_solidus
  S_cotti
  L_intestinalis
  Spirometra_Aus1
  D_latus
  S_proliferum
  E_multilocularis
  H_microstoma
  E_granulosus
  T_solium
  T_multiceps
)

FOCAL_SPECIES=(S_cotti S_solidus)

COTTI_RELEASE_ROOT=${COTTI_RELEASE_ROOT:-/work/shared/Scotti_release_v1.0_preNCBI_20261118}
COTTI_RELEASE_SHA256=${COTTI_RELEASE_SHA256:-07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371}
COTTI_PROTEINS=${COTTI_PROTEINS:-"${COTTI_RELEASE_ROOT}/proteins/Scotti_proteins.longest_isoform.v1.0.fa"}
# Raw/all-isoform input is distinct from the 9,234 representatives above.
COTTI_ALL_PROTEINS=${COTTI_ALL_PROTEINS:-"${COTTI_RELEASE_ROOT}/proteins/Scotti_proteins.all_isoforms.v1.0.fa"}
# Historical source root retained for scripts that explicitly request it.
COTTI_CANDIDATE_ROOT=${COTTI_CANDIDATE_ROOT:-/work/cyu/annotation/scotti-annotation/results/filtered_candidates/Scotti_filtered7_v1}
COTTI_GFF=${COTTI_GFF:-"${COTTI_RELEASE_ROOT}/annotation/Scotti_nuclear.protein_coding.v1.0.gff3"}
COTTI_GENOME=${COTTI_GENOME:-"${COTTI_RELEASE_ROOT}/assembly/Scotti_nuclear.fa"}

# Deliberately excluded from the primary analysis. PRJEB1202 is highly
# fragmented and its protein count is likely inflated; PRJEB35375 lacks a
# usable public structural annotation.
EXCLUDED_SPECIES=(S_erinaceieuropaei_PRJEB1202)

PENDING_SPECIES=()
