#!/usr/bin/env bash
# Copy the validated scientific release to /work/shared without overwriting an
# existing release. This is a pre-NCBI package, not an accessioned submission.
set -euo pipefail

ASSEMBLY=/work/cyu/assembly/worm/final_release
ANNOT=/work/cyu/annotation/ssol-annotation/results/final_filtered_release
FUNC=${ANNOT}/functional_integration
SHARED=${SHARED:-/work/shared/Ssolidus_release_v1.0_preNCBI_20260831}

[[ ! -e "${SHARED}" ]] || { echo "ERROR: destination exists: ${SHARED}" >&2; exit 1; }
mkdir -p "${SHARED}"/{assembly,annotation,proteins,functional_annotation,quality_control}

cp "${ASSEMBLY}/Ssolidus_nuclear.fa" "${ASSEMBLY}/Ssolidus_mitochondrial.fa" \
    "${ASSEMBLY}/Ssolidus_mitochondrial.gff3" "${SHARED}/assembly/"
cp "${ANNOT}/Ssolidus_nuclear.annotation.final_candidate.gff3" \
    "${SHARED}/annotation/Ssolidus_nuclear.annotation.v1.0.gff3"
cp "${ANNOT}/Ssolidus_nuclear.protein_coding.release_candidate.gtf" \
    "${SHARED}/annotation/Ssolidus_nuclear.protein_coding.v1.0.gtf"
cp "${ANNOT}/Ssolidus_nuclear.proteins.filtered62.fa" \
    "${SHARED}/proteins/Ssolidus_proteins.all_isoforms.v1.0.fa"
cp "${ANNOT}/Ssolidus_nuclear.filtered62.longest_isoform.fa" \
    "${SHARED}/proteins/Ssolidus_proteins.longest_isoform.v1.0.fa"
cp "${FUNC}/Ssolidus_protein_function_master.tsv" "${FUNC}/Ssolidus_product_name_assignment.tsv" \
    "${SHARED}/functional_annotation/"
cp "${ANNOT}/Ssolidus_nuclear.annotation.final_candidate.statistics.txt" \
    "${ANNOT}/release_candidate_coordinate_audit.final.tsv" \
    "${SHARED}/quality_control/"

(cd "${SHARED}" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
chgrp -R shared "${SHARED}" 2>/dev/null || true
chmod -R g+rwX "${SHARED}"
(cd "${SHARED}" && sha256sum -c SHA256SUMS)
echo "PASS: shared scientific release created at ${SHARED}"
