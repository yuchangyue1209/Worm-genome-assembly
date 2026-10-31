#!/usr/bin/env bash
# Merge functional protein-coding and conservative ncRNA annotations, normalize
# with AGAT, and remove whitespace inherited from tRNAscan-SE tabular fields.
set -euo pipefail

ASM=/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
FINAL=/work/cyu/annotation/ssol-annotation/results/final_filtered_release
PCGFF=${FINAL}/Ssolidus_nuclear.protein_coding.functional.gff3
NCRNA=/work/cyu/annotation/ssol-annotation/results/ncrna_final/release/Ssolidus_nuclear.ncRNA.v1.0.gff3
RAW=${FINAL}/Ssolidus_nuclear.annotation.functional_plus_ncRNA.raw.gff3
MERGED=${FINAL}/Ssolidus_nuclear.annotation.functional_plus_ncRNA.gff3
SANITIZED=${FINAL}/Ssolidus_nuclear.annotation.final_candidate.gff3

for cmd in samtools agat_convert_sp_gxf2gxf.pl; do
    command -v "${cmd}" >/dev/null || { echo "ERROR: ${cmd} not found" >&2; exit 1; }
done
for file in "${ASM}" "${PCGFF}" "${NCRNA}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

samtools faidx "${ASM}"
{
    echo '##gff-version 3'
    awk 'BEGIN{OFS=" "}{print "##sequence-region",$1,1,$2}' "${ASM}.fai"
    awk '/^##FASTA/{exit} $0!~/^#/&&NF{print}' "${PCGFF}"
    awk '/^##FASTA/{exit} $0!~/^#/&&NF{print}' "${NCRNA}"
} > "${RAW}"

agat_convert_sp_gxf2gxf.pl --gff "${RAW}" --output "${MERGED}"

awk -F'\t' 'BEGIN{OFS="\t"}
    /^#/{print;next}
    NF==9{for(i=1;i<=8;i++){sub(/^[[:space:]]+/,"",$i);sub(/[[:space:]]+$/, "",$i)};print;next}
    {print}
' "${MERGED}" > "${SANITIZED}"

genes=$(awk -F'\t' '$3=="gene"{n++}END{print n+0}' "${SANITIZED}")
transcripts=$(awk -F'\t' '$3=="transcript"{n++}END{print n+0}' "${SANITIZED}")
ncrna=$(awk -F'\t' '$3~/^(tRNA|rRNA|snRNA|snoRNA|miRNA|RNase_P_RNA)$/{n++}END{print n+0}' "${SANITIZED}")
[[ "${genes}" -eq 10945 && "${transcripts}" -eq 12849 && "${ncrna}" -eq 1149 ]] || {
    echo "ERROR: merged feature counts are incorrect" >&2; exit 1;
}
echo "PASS: merged annotation ${SANITIZED}"
