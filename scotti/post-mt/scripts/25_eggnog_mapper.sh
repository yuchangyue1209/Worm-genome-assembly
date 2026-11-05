#!/usr/bin/env bash
# Environment: eggnog_mapper. Corrected installed CLI: -i and -m, NOT --input/--method.
# First attempt stopped during argument parsing. Corrected execution result not yet reported.
# Accepts the existing directory with logs only; refuses existing emapper result files.
source "$(dirname "$0")/00_config.sh"
PROT="$ANNO/results/braker3_etp_v1_retry2/braker.aa"
DB=/work/cyu/databases/eggnog_mapper
OUT="$ANNO/results/functional/02_eggnog"
command -v emapper.py >/dev/null
for file in "$PROT" "$DB/eggnog.db" "$DB/eggnog.taxa.db" "$DB/eggnog_proteins.dmnd";do [[ -s "$file" ]];done
if compgen -G "$OUT/scotti_braker.emapper.*" > /dev/null;then
 echo 'Existing emapper outputs found. Inspect running processes/results before restarting.';exit 1
fi
mkdir -p "$OUT"
emapper.py --version > "$OUT/emapper.version.current.txt" 2>&1
emapper.py --help > "$OUT/emapper.help.current.txt" 2>&1
sha256sum "$PROT" > "$OUT/query.sha256"
LOG="$OUT/run.$(date +%Y%m%d_%H%M%S).log"
emapper.py --data_dir "$DB" -i "$PROT" -o scotti_braker --output_dir "$OUT" \
 --itype proteins -m diamond --cpu 32 2>&1 | tee "$LOG"
[[ -s "$OUT/scotti_braker.emapper.annotations" ]]
awk -F '\t' '!/^#/ && NF {rows++;ids[$1]=1} END {for(id in ids)n++;print "Annotation rows:",rows+0;print "Annotated proteins:",n+0}' "$OUT/scotti_braker.emapper.annotations"
# No --override and no solidus-specific fixed expected annotation count.
