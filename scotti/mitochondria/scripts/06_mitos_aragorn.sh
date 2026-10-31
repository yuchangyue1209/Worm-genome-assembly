#!/usr/bin/env bash
# Environment: mitos2. Real entry point is runmitos (not runmitos.py).
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
need "$ASM"
command -v runmitos
[[ -d "$DBBASE/refseq89m" ]] || { echo "Missing database: $DBBASE/refseq89m"; exit 1; }
fresh "$MITOS"
runmitos -i "$ASM" -o "$MITOS" -R "$DBBASE" -r refseq89m -c 9 --best 2>&1 | tee "$MITOS/run.log"
/home/cyu/.conda/envs/genome_assembly/bin/aragorn -mt -gc9 "$ASM" > "$MITOS/aragorn.scotti.out"
# Print results; preserve all original predictions, including OH/OL.
find "$MITOS" -maxdepth 2 -type f -name '*.gff*' -exec awk -F '\t' '/^##FASTA/{exit} /^#/{next} NF>=9 {print $3 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' {} \;
grep -n -F -B 40 -A 4 'mtRNA-His' "$MITOS/aragorn.scotti.out" || true
