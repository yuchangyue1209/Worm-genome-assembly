#!/usr/bin/env bash
# Validate the current FCS-clean assembly and completed annotation components.
# Functional annotation is not final until eggNOG/InterPro/product naming and
# NCBI validation have been completed.
set -euo pipefail

ASM=/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
MT=/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa
MTGFF=/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.gff3
GFF=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3
GTF=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.checked.gtf
PROT=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_final.proteins.fa
NCRNA=/work/cyu/annotation/ssol-annotation/results/ncrna_final/release/Ssolidus_nuclear.ncRNA.v1.0.gff3
SWISS=/work/cyu/annotation/ssol-annotation/results/functional_final/swissprot/ssol_final_vs_swissprot.best_hit.tsv
OUT=/work/cyu/annotation/ssol-annotation/results/release_qc/final_manifest
mkdir -p "${OUT}"

FILES=("${ASM}" "${MT}" "${MTGFF}" "${GFF}" "${GTF}" "${PROT}" "${NCRNA}" "${SWISS}")
for file in "${FILES[@]}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

sha256sum "${FILES[@]}" > "${OUT}/SHA256SUMS.current"

seqs=$(grep -c '^>' "${ASM}")
bases=$(awk '/^>/{next}{gsub(/[[:space:]]/,"");n+=length}END{print n+0}' "${ASM}")
genes=$(awk -F'\t' '$3=="gene"{n++}END{print n+0}' "${GFF}")
transcripts=$(awk -F'\t' '$3=="mRNA"||$3=="transcript"{n++}END{print n+0}' "${GFF}")
proteins=$(grep -c '^>' "${PROT}")
ncrnas=$(awk -F'\t' '$0!~/^#/{n++}END{print n+0}' "${NCRNA}")
swiss=$(awk -F'\t' 'NF&&!seen[$1]++{n++}END{print n+0}' "${SWISS}")

cat > "${OUT}/release_status.tsv" <<EOF
metric	value
nuclear_assembly_bp	${bases}
nuclear_scaffolds	${seqs}
protein_coding_genes	${genes}
transcripts_or_isoforms	${transcripts}
proteins	${proteins}
conservative_nuclear_ncRNA	${ncrnas}
proteins_with_SwissProt_hit	${swiss}
functional_release_status	pending_eggNOG_InterPro_product_names_and_NCBI_validation
EOF

[[ "${bases}" -eq 846268314 && "${seqs}" -eq 166 ]] || {
    echo "ERROR: assembly does not match the FCS-clean release" >&2; exit 1;
}
[[ "${genes}" -eq 9858 && "${proteins}" -eq 12911 && "${ncrnas}" -eq 1149 ]] || {
    echo "ERROR: annotation counts differ from the documented release" >&2; exit 1;
}

column -t -s $'\t' "${OUT}/release_status.tsv"
echo "Checksums: ${OUT}/SHA256SUMS.current"
