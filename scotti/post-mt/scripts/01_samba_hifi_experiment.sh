#!/usr/bin/env bash
# Environment: genome_assembly. Historical experiment; NOT selected annotation input.
source "$(dirname "$0")/00_config.sh"
OUT="${OUT:-$ROOT/scaffolding/samba_hifi_asm_m5000_reconstructed}"
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"; cd "$OUT"
python3 - "$HIFI" "$OUT/hifi.fa" <<'EOF'
import gzip, sys
from Bio import SeqIO
n=b=0
with gzip.open(sys.argv[1],'rt') as src, open(sys.argv[2],'x') as dst:
    for r in SeqIO.parse(src,'fastq'):
        dst.write('>'+r.id+'\n'+str(r.seq)+'\n'); n+=1; b+=len(r)
if (n,b)!=(4662246,26074063243):
    raise SystemExit(f'Unexpected input: {n} reads, {b} bases')
print(n,b)
EOF
cp "$RAGTAG" input.fa
samba.sh -r input.fa -q "$OUT/hifi.fa" -d asm -t 30 -m 5000 2>&1 | tee run.log
