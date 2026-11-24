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

# Ten-species background set. Keep stable standardized inputs unchanged.
# Keep this background array stable for existing functional/secretome workflows.
# cotti release inclusion is explicit in script 03; it does not silently alter other jobs.
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

# cotti is released; standardized links in an existing server project remain a separate action.
PENDING_SPECIES=()
RELEASED_FOCAL_SPECIES=(S_solidus S_cotti)
SCOTTI_RELEASE=${SCOTTI_RELEASE:-/work/shared/Scotti_release_v1.0_preNCBI_20261118}
INCLUDE_SCOTTI_RELEASE=${INCLUDE_SCOTTI_RELEASE:-0}

# Six-taxon core describes the biological design; it is not a command to remove outgroups.
FOCAL_SPECIES=(S_solidus S_cotti)
CORE_DIPHYLLOBOTHRIIDEA=(S_solidus S_cotti L_intestinalis D_latus Spirometra_Aus1 S_proliferum)
# Historical candidate mode is retained for reproducing old input freezes.
INCLUDE_SCOTTI_CANDIDATE=${INCLUDE_SCOTTI_CANDIDATE:-0}
SCOTTI_CANDIDATE=${SCOTTI_CANDIDATE:-/work/cyu/annotation/scotti-annotation/results/filtered_candidates/Scotti_filtered7_v1}
SCOTTI_CANDIDATE_STAGE=${SCOTTI_CANDIDATE_STAGE:-"${PROJECT}/02_standardized_candidates/Scotti_filtered7_v1"}
# Unique new run ID prevents stale symlinks/taxa from an older OrthoFinder input directory.
ORTHO_RUN_ID=${ORTHO_RUN_ID:-"background10_$(date +%Y%m%d_%H%M%S)"}
