#!/usr/bin/env bash

SCOTTI_ROOT=${SCOTTI_ROOT:-/work/cyu/annotation/scotti-annotation}
BRAKER_DIR=${BRAKER_DIR:-"${SCOTTI_ROOT}/results/braker3_etp_v1_retry2"}
PROTEINS=${PROTEINS:-"${BRAKER_DIR}/braker.aa"}
GTF=${GTF:-"${BRAKER_DIR}/braker.gtf"}
CDS=${CDS:-"${BRAKER_DIR}/braker.codingseq"}

FUNCTIONAL=${FUNCTIONAL:-"${SCOTTI_ROOT}/results/functional"}
EGGNOG_OUT=${EGGNOG_OUT:-"${FUNCTIONAL}/02_eggnog"}
INTERPRO_OUT=${INTERPRO_OUT:-"${FUNCTIONAL}/03_interproscan"}

EGGNOG_DB=${EGGNOG_DB:-/work/cyu/databases/eggnog_mapper}
INTERPROSCAN=${INTERPROSCAN:-/work/cyu/software/interproscan-5.78-109.0/interproscan.sh}
THREADS=${THREADS:-32}
EXPECTED_PROTEINS=${EXPECTED_PROTEINS:-11471}
