#!/usr/bin/env bash
# Environment: genome_assembly. Validation completed historically; inspect existing logs first.
source "$(dirname "$0")/00_config.sh"
DIR="$ANNO/results/filtered_candidates/Scotti_filtered7_v1"
GT=/home/cyu/.conda/envs/genome_assembly/bin/gt
# Script 27 now writes ##gff-version 3 first. The original server package needed a header-only fix.
"$GT" gff3validator -typecheck so "$DIR/Scotti_nuclear.filtered7.candidate.gff3" \
 2>&1 | tee "$DIR/genometools_validation.log"
if [[ -e "$DIR/busco_longest_odb12" ]];then
 echo 'Existing BUSCO run retained; review its summary instead of rerunning.'
else
 busco -i "$DIR/Scotti_nuclear.filtered7.longest.proteins.fa" -m proteins \
  -l metazoa_odb12 --cpu 20 --out_path "$DIR" -o busco_longest_odb12 \
  2>&1 | tee "$DIR/busco_longest_odb12.log"
fi
# gff3validator warnings about absent sequence-region directives are not fatal;
# the validator does not rewrite the source file or validate translations against the genome.
