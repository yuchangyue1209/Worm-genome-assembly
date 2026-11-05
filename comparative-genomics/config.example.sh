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
BUSCO_LINEAGE=${BUSCO_LINEAGE:-metazoa_odb10}
THREADS=${THREADS:-32}

# Current analysis set. The S_cotti BRAKER candidate and protein BUSCO runs are
# complete, but it must be added only after final filtering, release and the
# same standardization checks.
SPECIES=(
  S_solidus
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

# Deliberately excluded from the primary analysis. PRJEB1202 is highly
# fragmented and its protein count is likely inflated; PRJEB35375 lacks a
# usable public structural annotation.
EXCLUDED_SPECIES=(S_erinaceieuropaei_PRJEB1202)

# Add to SPECIES only when curated genome, GFF3, CDS, and primary proteins are
# released and linked into the standardized project layout.
PENDING_SPECIES=(S_cotti)
