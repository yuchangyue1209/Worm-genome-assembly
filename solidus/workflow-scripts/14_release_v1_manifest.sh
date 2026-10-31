#!/usr/bin/env bash
# Validate the completed filtered structural, functional, and ncRNA release.
set -euo pipefail

ASM=/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
MT=/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa
MTGFF=/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.gff3
ROOT=/work/cyu/annotation/ssol-annotation/results/final_filtered_release
GFF=${ROOT}/Ssolidus_nuclear.annotation.final_candidate.gff3
GTF=${ROOT}/Ssolidus_nuclear.protein_coding.release_candidate.gtf
PROT=${ROOT}/Ssolidus_nuclear.proteins.filtered62.fa
LONGEST=${ROOT}/Ssolidus_nuclear.filtered62.longest_isoform.fa
MASTER=${ROOT}/functional_integration/Ssolidus_protein_function_master.tsv
PRODUCTS=${ROOT}/functional_integration/Ssolidus_product_name_assignment.tsv
OUT=${ROOT}/release_manifest
mkdir -p "${OUT}"

FILES=("${ASM}" "${MT}" "${MTGFF}" "${GFF}" "${GTF}" "${PROT}" "${LONGEST}" "${MASTER}" "${PRODUCTS}")
for file in "${FILES[@]}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

seqs=$(grep -c '^>' "${ASM}")
bases=$(awk '/^>/{next}{gsub(/[[:space:]]/,"");n+=length}END{print n+0}' "${ASM}")
pc_genes=$(awk -F'\t' '$3=="gene"{n++}END{print n+0}' \
    /work/cyu/annotation/ssol-annotation/results/final_filtered_release/Ssolidus_nuclear.protein_coding.functional.gff3)
genes=$(awk -F'\t' '$3=="gene"{n++}END{print n+0}' "${GFF}")
transcripts=$(awk -F'\t' '$3=="transcript"||$3=="mRNA"{n++}END{print n+0}' "${GFF}")
proteins=$(grep -c '^>' "${PROT}")
ncrna=$((genes-pc_genes))
assigned=$(awk -F'\t' 'NR>1&&$4!="hypothetical protein"{n++}END{print n+0}' "${PRODUCTS}")
hypothetical=$(awk -F'\t' 'NR>1&&$4=="hypothetical protein"{n++}END{print n+0}' "${PRODUCTS}")

printf '%s\t%s\n' \
    metric value \
    nuclear_assembly_bp "${bases}" \
    nuclear_scaffolds "${seqs}" \
    protein_coding_genes "${pc_genes}" \
    protein_coding_transcripts "${transcripts}" \
    protein_isoforms "${proteins}" \
    conservative_nuclear_ncRNA "${ncrna}" \
    total_gene_features "${genes}" \
    proteins_with_assigned_product "${assigned}" \
    hypothetical_proteins "${hypothetical}" \
    scientific_release_status complete_preNCBI \
    NCBI_submission_status pending_locus_tag_prefix_and_validator \
    > "${OUT}/release_status.tsv"

[[ "${bases}" -eq 846268314 && "${seqs}" -eq 166 ]] || {
    echo "ERROR: assembly differs from the FCS-clean release" >&2; exit 1;
}
[[ "${pc_genes}" -eq 9796 && "${transcripts}" -eq 12849 && "${proteins}" -eq 12849 ]] || {
    echo "ERROR: protein-coding counts differ from the filtered release" >&2; exit 1;
}
[[ "${ncrna}" -eq 1149 && "${genes}" -eq 10945 ]] || {
    echo "ERROR: ncRNA or total-gene count differs from the final release" >&2; exit 1;
}
[[ "${assigned}" -eq 10535 && "${hypothetical}" -eq 2314 ]] || {
    echo "ERROR: product-name counts differ from the reviewed release" >&2; exit 1;
}

sha256sum "${FILES[@]}" > "${OUT}/SHA256SUMS"
column -t -s $'\t' "${OUT}/release_status.tsv"
echo "PASS: completed scientific release manifest"
