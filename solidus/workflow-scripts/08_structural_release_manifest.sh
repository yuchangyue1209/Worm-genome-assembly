#!/usr/bin/env bash
# Validate the final nuclear/mitochondrial assembly and annotation release.
set -euo pipefail

RELEASE=/work/cyu/assembly/worm/final_release
ANNOTATION=/work/cyu/annotation/ssol-annotation/results/tsa_alignment

FILES=(
    "${RELEASE}/Ssolidus_nuclear.fa"
    "${RELEASE}/Ssolidus_mitochondrial.fa"
    "${RELEASE}/Ssolidus_mitochondrial.gff3"
    "${ANNOTATION}/ssol_braker_tsa_utr_tagseq_plus_187.checked.gtf"
    "${ANNOTATION}/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3"
)

for file in "${FILES[@]}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done

echo "Release files"
for file in "${FILES[@]}"; do
    printf '%s\t' "${file}"
    sha256sum "${file}" | cut -d' ' -f1
done

awk '
    /^>/ {
        if (seen) {n++; total+=len; if(len>max) max=len; lengths[n]=len}
        seen=1; len=0; next
    }
    {gsub(/[[:space:]]/, ""); len+=length($0)}
    END {
        if (seen) {n++; total+=len; if(len>max) max=len; lengths[n]=len}
        for(i=1;i<=n;i++) for(j=i+1;j<=n;j++) if(lengths[j]>lengths[i]) {
            t=lengths[i]; lengths[i]=lengths[j]; lengths[j]=t
        }
        cumulative=0
        for(i=1;i<=n;i++) {cumulative+=lengths[i]; if(cumulative>=total/2) {n50=lengths[i]; break}}
        printf "nuclear_bp\t%d\nscaffolds\t%d\nN50_bp\t%d\nmax_bp\t%d\n", total,n,n50,max
    }
' "${RELEASE}/Ssolidus_nuclear.fa"

gff="${ANNOTATION}/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3"
for feature in gene mRNA exon CDS five_prime_UTR three_prime_UTR; do
    awk -F '\t' -v feature="${feature}" '$3==feature{n++} END{print feature "\t" n+0}' "${gff}"
done
