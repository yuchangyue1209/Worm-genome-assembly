#!/usr/bin/env bash
# PENDING template adapted from solidus. Not yet run or validated against installed CLI.
# Environment must provide the Java version required by this InterProScan installation.
# Do not launch until resources and command help have been reviewed.
source "$(dirname "$0")/00_config.sh"
[[ "${RUN_INTERPRO:-0}" == 1 ]] || { echo 'Pending template. Review before setting RUN_INTERPRO=1.';exit 0; }
IPR=/work/cyu/software/interproscan-5.78-109.0
OUT="$ANNO/results/functional/03_interproscan"
PROT="$ANNO/results/braker3_etp_v1_retry2/braker.aa"
[[ -x "$IPR/interproscan.sh" && -s "$PROT" ]]
command -v java
java -version
mkdir "$OUT"
python3 - "$PROT" "$OUT/scotti_braker.no_terminal_stop.fa" <<'EOF'
import sys
from Bio import SeqIO
allowed=set('ACDEFGHIKLMNPQRSTVWYX')
records=[]
for r in SeqIO.parse(sys.argv[1],'fasta'):
    seq=str(r.seq).upper()
    if seq.endswith('*'):seq=seq[:-1]
    if not seq or set(seq)-allowed:raise ValueError(f'Unexpected residues in {r.id}: {set(seq)-allowed}')
    r.seq=type(r.seq)(seq);r.description='';records.append(r)
with open(sys.argv[2],'x') as f:SeqIO.write(records,f,'fasta')
print('Prepared proteins:',len(records))
EOF
# Strip only a terminal stop; preserve IDs and X, reject unsupported residues.
"$IPR/interproscan.sh" -i "$OUT/scotti_braker.no_terminal_stop.fa" \
 -b "$OUT/scotti_braker" -f TSV,GFF3 -goterms -pa -cpu 32 \
 2>&1 | tee "$OUT/run.log"
