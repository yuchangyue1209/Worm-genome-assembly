#!/usr/bin/env bash
# Build the conservative nuclear ncRNA release from completed Rfam/tRNAscan
# results and audit overlap with protein-coding CDS. This script does not rerun
# the four-hour whole-genome cmscan.
set -euo pipefail

ROOT=/work/cyu/annotation/ssol-annotation/results/ncrna_final
RFAM=${ROOT}/rfam_genome
TRNA_IN=${ROOT}/trna/complete/Ssolidus_tRNAscan.out
RFGFF=${RFAM}/Ssolidus_nuclear.rfam.conservative.v2.gff3
TRNA_DIR=${ROOT}/trna/conservative_functional
TRNA_GFF=${TRNA_DIR}/Ssolidus_tRNA.conservative.gff3
FINAL=${ROOT}/Ssolidus_nuclear.ncRNA.conservative.gff3
PCGFF=/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3
AUDIT=${ROOT}/overlap_audit
RELEASE=${ROOT}/release

for cmd in bedtools awk sort sha256sum; do
    command -v "${cmd}" >/dev/null || { echo "ERROR: ${cmd} not found" >&2; exit 1; }
done
for file in "${TRNA_IN}" "${RFGFF}" "${PCGFF}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done
mkdir -p "${TRNA_DIR}" "${AUDIT}" "${RELEASE}"

awk -F'\t' '
BEGIN {OFS="\t"; print "##gff-version 3"}
NR>3 && NF>=11 && tolower($12) !~ /pseudo/ && $5==$10 {
    begin=$3+0; end=$4+0
    if(begin<=end){start=begin; stop=end; strand="+"}
    else{start=end; stop=begin; strand="-"}
    id=sprintf("ssol_tRNA_%04d",++n)
    print $1,"tRNAscan-SE","tRNA",start,stop,$9+0,strand,".",
      "ID=" id ";Name=tRNA-" $5 "-" $6 ";product=tRNA-" $5 \
      ";anticodon=" $6 ";isotype=" $5 ";domain_score=" ($9+0) \
      ";isotype_score=" ($11+0) ";tRNAscan_number=" $2
}' "${TRNA_IN}" > "${TRNA_GFF}"

trna_n=$(awk -F'\t' '$0!~/^#/&&NF>=9{n++}END{print n+0}' "${TRNA_GFF}")
[[ "${trna_n}" -eq 352 ]] || { echo "ERROR: expected 352 tRNAs, found ${trna_n}" >&2; exit 1; }

{
    echo '##gff-version 3'
    grep -v '^#' "${RFGFF}"
    grep -v '^#' "${TRNA_GFF}"
} | sort -t$'\t' -k1,1V -k4,4n -k5,5n > "${FINAL}"

total=$(awk -F'\t' '$0!~/^#/&&NF>=9{n++}END{print n+0}' "${FINAL}")
[[ "${total}" -eq 1149 ]] || { echo "ERROR: expected 1,149 ncRNAs, found ${total}" >&2; exit 1; }

awk -F'\t' 'BEGIN{OFS="\t"}$0!~/^#/&&$3=="CDS"{print $1,($4+0)-1,$5+0,$9,".",$7}' \
    "${PCGFF}" | sort -k1,1V -k2,2n > "${AUDIT}/protein_coding_CDS.bed"
awk -F'\t' 'BEGIN{OFS="\t"}$0!~/^#/&&NF>=9{print $1,($4+0)-1,$5+0,$3"|"$9,$6,$7}' \
    "${FINAL}" | sort -k1,1V -k2,2n > "${AUDIT}/conservative_ncRNA.bed"

bedtools intersect -a "${AUDIT}/conservative_ncRNA.bed" \
    -b "${AUDIT}/protein_coding_CDS.bed" -wao > "${AUDIT}/ncRNA_vs_CDS.all_overlaps.tsv"
bedtools intersect -s -a "${AUDIT}/conservative_ncRNA.bed" \
    -b "${AUDIT}/protein_coding_CDS.bed" -wao > "${AUDIT}/ncRNA_vs_CDS.same_strand.tsv"

all=$(awk -F'\t' '$13>0{n++}END{print n+0}' "${AUDIT}/ncRNA_vs_CDS.all_overlaps.tsv")
same=$(awk -F'\t' '$13>0{n++}END{print n+0}' "${AUDIT}/ncRNA_vs_CDS.same_strand.tsv")
cat > "${AUDIT}/ncRNA_CDS_overlap_conclusion.tsv" <<EOF
metric	value
conservative_ncRNA_features	${total}
all_strands_CDS_overlap_records	${all}
same_strand_CDS_overlap_records	${same}
decision	retain_all_conservative_ncRNA_features
EOF

cp "${FINAL}" "${RELEASE}/Ssolidus_nuclear.ncRNA.v1.0.gff3"
cp "${RFAM}/Ssolidus_nuclear.rfam.comprehensive.gff3" \
    "${RELEASE}/Ssolidus_nuclear.ncRNA.comprehensive_candidates.v1.0.gff3"
sha256sum "${RELEASE}"/*.gff3 > "${RELEASE}/SHA256SUMS"

awk -F'\t' '$0!~/^#/{count[$3]++;total++}END{for(x in count)print x,count[x];print "TOTAL",total}' \
    "${FINAL}" | sort | column -t
column -t -s $'\t' "${AUDIT}/ncRNA_CDS_overlap_conclusion.tsv"
