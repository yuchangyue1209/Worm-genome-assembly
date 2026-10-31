#!/usr/bin/env bash
# Final sequence-ID, coordinate, count, AGAT, and checksum validation.
set -euo pipefail

ASM=/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
FINAL=/work/cyu/annotation/ssol-annotation/results/final_filtered_release
GFF=${FINAL}/Ssolidus_nuclear.annotation.final_candidate.gff3
GTF=${FINAL}/Ssolidus_nuclear.protein_coding.release_candidate.gtf
PROT=${FINAL}/Ssolidus_nuclear.proteins.filtered62.fa
LONGEST=${FINAL}/Ssolidus_nuclear.filtered62.longest_isoform.fa
MASTER=${FINAL}/functional_integration/Ssolidus_protein_function_master.tsv
PRODUCTS=${FINAL}/functional_integration/Ssolidus_product_name_assignment.tsv
AUDIT=${FINAL}/release_candidate_coordinate_audit.final.tsv

for file in "${ASM}" "${GFF}" "${GTF}" "${PROT}" "${LONGEST}" "${MASTER}" "${PRODUCTS}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done
samtools faidx "${ASM}"

awk -F'\t' 'BEGIN{OFS="\t";print "problem","sequence","feature","start","end","detail"}
NR==FNR{len[$1]=$2;next}/^#/{next}
NF!=9{print "invalid_column_count",$1,$3,$4,$5,NF;next}
!($1 in len){print "unknown_sequence",$1,$3,$4,$5,"not_in_assembly"}
$4!~/^[0-9]+$/||$5!~/^[0-9]+$/{print "non_integer_coordinate",$1,$3,$4,$5,"invalid";next}
$4<1{print "start_less_than_1",$1,$3,$4,$5,"invalid"}
$5<$4{print "end_before_start",$1,$3,$4,$5,"invalid"}
($1 in len)&&$5>len[$1]{print "end_beyond_sequence",$1,$3,$4,$5,len[$1]}
$7!~/^[-+.?]$/{print "invalid_strand",$1,$3,$4,$5,$7}
$3=="CDS"&&$8!~/^[012]$/{print "invalid_CDS_phase",$1,$3,$4,$5,$8}
$3!="CDS"&&$8!="."{print "unexpected_nonCDS_phase",$1,$3,$4,$5,$8}
' "${ASM}.fai" "${GFF}" > "${AUDIT}"

[[ "$(wc -l < "${AUDIT}")" -eq 1 ]] || { column -t -s $'\t' "${AUDIT}" | head -30; exit 1; }
agat_sp_statistics.pl --gff "${GFF}" --output "${FINAL}/Ssolidus_nuclear.annotation.final_candidate.statistics.txt"
sha256sum "${ASM}" "${GFF}" "${GTF}" "${PROT}" "${LONGEST}" "${MASTER}" "${PRODUCTS}" \
    > "${FINAL}/Ssolidus_final_candidate.SHA256SUMS"
echo "PASS: final scientific release validation completed"
