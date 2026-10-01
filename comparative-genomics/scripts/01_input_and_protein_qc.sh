#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
require_command seqkit

OUT=${QC}/comparative_input_qc
mkdir -p "${OUT}"

for sp in "${SPECIES[@]}"; do
    require_file "${STD}/genomes/${sp}.fasta"
    require_file "${STD}/proteins_primary/${sp}.faa"
done

genomes=()
proteins=()
for sp in "${SPECIES[@]}"; do
    genomes+=("${STD}/genomes/${sp}.fasta")
    proteins+=("${STD}/proteins_primary/${sp}.faa")
done

seqkit stats -a -T "${genomes[@]}" > "${OUT}/assembly_seqkit_stats.tsv"
seqkit stats -a -T "${proteins[@]}" > "${OUT}/protein_seqkit_stats.tsv"

printf 'Species\tProteins\tMean_length\tMedian_length\tMin_length\tMax_length\tLT30\tLT50\tLT100\tGT2000\tDuplicate_IDs\n' \
    > "${OUT}/protein_length_qc.tsv"

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    lengths=${OUT}/${sp}.lengths.txt
    awk '
        /^>/ {if (seen) print length(seq); seq=""; seen=1; next}
        {gsub(/[[:space:]]/, ""); seq=seq $0}
        END {if (seen) print length(seq)}
    ' "${protein}" | sort -n > "${lengths}"

    duplicate_ids=$(awk '/^>/{id=$1; sub(/^>/,"",id); seen[id]++} END{for(id in seen) if(seen[id]>1)n++; print n+0}' "${protein}")
    awk -v sp="${sp}" -v duplicate_ids="${duplicate_ids}" '
        {a[NR]=$1; sum+=$1; if($1<30)x30++; if($1<50)x50++; if($1<100)x100++; if($1>2000)x2000++}
        END {
            n=NR; med=(n%2 ? a[(n+1)/2] : (a[n/2]+a[n/2+1])/2)
            printf "%s\t%d\t%.1f\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", sp,n,sum/n,med,a[1],a[n],x30+0,x50+0,x100+0,x2000+0,duplicate_ids
        }
    ' "${lengths}" >> "${OUT}/protein_length_qc.tsv"
done

echo "PASS: input and protein QC written to ${OUT}"
