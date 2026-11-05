#!/usr/bin/env bash
# Environment: rnaseq_annotation. STAR 2.7.11b. Reconstructs the successful index.
source "$(dirname "$0")/00_config.sh"
BASE="$ANNO/results/star"
REFDIR="$BASE/01_index"
OUT="$BASE/01_index_no_gtf"
mkdir -p "$REFDIR"
REF="$REFDIR/reference.nuclear_plus_mt.fa"
# Preserve the reference left by the historical failed index attempt if it exists.
if [[ ! -e "$REF" ]]; then
python3 - "$NUCLEAR" "$MT" "$REF" <<'EOF'
from pathlib import Path
import sys
nuclear,mt,out=map(Path,sys.argv[1:]); ids=set(); totals=[0,0]
with out.open('x') as dst:
    for path,is_mt in ((nuclear,False),(mt,True)):
        with path.open() as src:
            for line in src:
                if line.startswith('>'):
                    name='Scotti_mt' if is_mt else line[1:].split()[0]
                    if name in ids: raise ValueError(f'Duplicate reference ID: {name}')
                    ids.add(name); totals[0]+=1; dst.write('>'+name+'\n')
                else:
                    totals[1]+=len(line.strip()); dst.write(line)
assert totals==[413,991332483], totals
EOF
fi
mkdir "$OUT"
# No --sjdbOverhang when building this annotation-free index: this version rejects it.
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir "$OUT" \
  --genomeFastaFiles "$REF" --outFileNamePrefix "$OUT/"
