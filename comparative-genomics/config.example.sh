#!/usr/bin/env bash

# Copy to config.sh and edit for the server. Do not commit private paths or data.
PROJECT=${PROJECT:-"${HOME}/Schistocephalus_comparative_genomics"}
RAW=${RAW:-"${PROJECT}/01_raw_data"}
STD=${STD:-"${PROJECT}/02_standardized_data"}
QC=${QC:-"${PROJECT}/03_qc"}
ORTHO_ROOT=${ORTHO_ROOT:-"${PROJECT}/04_orthofinder"}
FUNC=${FUNC:-"${PROJECT}/05_functional_annotation"}
FAMILY_QC=${FAMILY_QC:-"${PROJECT}/05_gene_family_qc"}
SECRETOME=${SECRETOME:-"${PROJECT}/06_secretome"}

INTERPROSCAN=${INTERPROSCAN:-/work/cyu/software/interproscan-5.78-109.0/interproscan.sh}
EGGNOG_DB=${EGGNOG_DB:-/work/cyu/databases/eggnog_mapper}
BUSCO_LINEAGE=${BUSCO_LINEAGE:-metazoa_odb10}
THREADS=${THREADS:-32}

# Current analysis set. S_cotti must be added only after its final annotation
# has passed the same standardization and QC steps.
SPECIES=(
  S_solidus
  L_intestinalis
  Spirometra_Aus1
  D_latus
  E_multilocularis
  H_microstoma
  E_granulosus
  T_solium
  T_multiceps
)

# Optional sensitivity-analysis taxon. PRJEB1202 is highly fragmented and
# should not be used for strict gene-loss or CAFE5 inference.
OPTIONAL_SPECIES=(S_erinaceieuropaei_PRJEB1202)

# Add to SPECIES only when genome, GFF3, CDS, and primary proteins are final.
PENDING_SPECIES=(S_cotti)
