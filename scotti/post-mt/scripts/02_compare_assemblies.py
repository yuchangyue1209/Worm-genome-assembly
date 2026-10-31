#!/usr/bin/env python3
"""Biopython required. Exact sequence comparison ignores names/order/orientation.
Usage: python 02_compare_assemblies.py OUTPUT_DIR label=/absolute/assembly.fa ...
Writes changed FASTAs for the first two assemblies, suitable for BLAST/PAF review.
"""
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import hashlib, itertools, sys
out=Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=False)
def digest(s):
    s=s.upper(); r=str(Seq(s).reverse_complement())
    return hashlib.sha256(min(s,r).encode()).hexdigest()
data={}; rows=['assembly\tsequences\ttotal_bp\tN_bp\tN50_bp\tlongest_bp']
for arg in sys.argv[2:]:
    label,path=arg.split('=',1)
    if not label.replace('_','').isalnum(): raise ValueError('Unsafe label')
    records=list(SeqIO.parse(path,'fasta')); lengths=sorted((len(r) for r in records),reverse=True)
    total=sum(lengths); cumulative=0; n50=0
    for size in lengths:
        cumulative+=size
        if cumulative*2>=total: n50=size; break
    counts=Counter(digest(str(r.seq)) for r in records)
    data[label]=(path,counts)
    rows.append(f'{label}\t{len(records)}\t{total}\t{sum(str(r.seq).upper().count("N") for r in records)}\t{n50}\t{max(lengths,default=0)}')
    del records
for a,b in itertools.combinations(data,2):
    ca,cb=data[a][1],data[b][1]
    rows.append(f'{a} vs {b}: shared={sum((ca&cb).values())}; identical_composition={ca==cb}')
labels=list(data)
if len(labels)>=2:
    for label,other in ((labels[0],labels[1]),(labels[1],labels[0])):
        remaining=data[other][1].copy()
        with open(out/f'{label}.changed.fa','x') as dst:
            for r in SeqIO.parse(data[label][0],'fasta'):
                key=digest(str(r.seq))
                if remaining[key]: remaining[key]-=1
                else: SeqIO.write(r,dst,'fasta')
(out/'assembly_comparison.txt').write_text('\n'.join(rows)+'\n')
print('\n'.join(rows))
