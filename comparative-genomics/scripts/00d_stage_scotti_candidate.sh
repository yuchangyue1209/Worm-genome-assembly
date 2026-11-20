#!/usr/bin/env bash
# Stages a separately labeled candidate. Does not replace STD/S_cotti or run OrthoFinder.
set -euo pipefail
source "$(dirname "$0")/common.sh"
CANDIDATE=${SCOTTI_CANDIDATE:-/work/cyu/annotation/scotti-annotation/results/filtered_candidates/Scotti_filtered7_v1}
STAGE=${SCOTTI_CANDIDATE_STAGE:-"${PROJECT}/02_standardized_candidates/Scotti_filtered7_v1"}
GENOME=/work/cyu/assembly/worm/scotti/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa
for file in "$GENOME" "$CANDIDATE/Scotti_nuclear.filtered7.candidate.gff3" \
 "$CANDIDATE/Scotti_nuclear.filtered7.cds.fa" "$CANDIDATE/Scotti_nuclear.filtered7.proteins.fa" \
 "$CANDIDATE/Scotti_nuclear.filtered7.longest.proteins.fa" "$CANDIDATE/summary.json";do require_file "$file";done
[[ ! -e "$STAGE" ]] || { echo 'Candidate stage exists; inspect rather than overwrite';exit 1; }
python3 - "$CANDIDATE" <<'EOF'
from pathlib import Path
import json,sys
p=Path(sys.argv[1]);s=json.loads((p/'summary.json').read_text())
assert (s['genes'],s['proteins'],s['longest_proteins'])==(9234,11464,9234)
assert s['status']=='candidate_not_final_release'
EOF
mkdir -p "$STAGE"/{genomes,gff,cds,proteins_all,proteins_primary}
ln -s "$GENOME" "$STAGE/genomes/S_cotti.fasta"
ln -s "$CANDIDATE/Scotti_nuclear.filtered7.candidate.gff3" "$STAGE/gff/S_cotti.gff3"
ln -s "$CANDIDATE/Scotti_nuclear.filtered7.cds.fa" "$STAGE/cds/S_cotti.fna"
ln -s "$CANDIDATE/Scotti_nuclear.filtered7.proteins.fa" "$STAGE/proteins_all/S_cotti.faa"
ln -s "$CANDIDATE/Scotti_nuclear.filtered7.longest.proteins.fa" "$STAGE/proteins_primary/S_cotti.faa"
printf '%s\n' 'EXPLORATORY CANDIDATE ONLY. Not a final filtered release. Stable standardized data were not replaced.' > "$STAGE/STATUS.txt"
sha256sum "$STAGE"/genomes/* "$STAGE"/gff/* "$STAGE"/cds/* "$STAGE"/proteins_all/* "$STAGE"/proteins_primary/* > "$STAGE/inputs.sha256"
echo "Candidate staged: $STAGE"
