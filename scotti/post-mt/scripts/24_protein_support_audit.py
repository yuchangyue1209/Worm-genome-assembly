#!/usr/bin/env python3
"""Reconstructed completed Swiss-Prot support audit. No automatic protein removal."""
from pathlib import Path
from Bio import SeqIO
from collections import Counter
import argparse,csv
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
root=Path('/work/cyu/annotation/scotti-annotation');a.out.mkdir(parents=True,exist_ok=False)
hits={}
for line in (root/'results/functional/01_swissprot/scotti_vs_swissprot.best_hit.tsv').open():
    f=line.rstrip('\n').split('	',9);hits[f[0]]=f
counts=Counter()
with (a.out/'protein_support.tsv').open('x') as f:
    w=csv.writer(f,delimiter='	');w.writerow(['protein_id','length_aa','contains_X','internal_stop','swissprot_hit','identity_pct','query_coverage_pct','subject_coverage_hsp_pct','evalue','bitscore','title'])
    for r in SeqIO.parse(root/'results/braker3_etp_v1_retry2/braker.aa','fasta'):
        seq=str(r.seq).upper();body=seq[:-1] if seq.endswith('*') else seq;n=len(body);h=hits.get(r.id)
        for label,flag in [('lt50aa',n<50),('gt5000aa',n>5000),('contains_X','X' in body)]:
            if flag:counts[label+'_total']+=1;counts[label+'_with_hit']+=h is not None
        w.writerow([r.id,n,'X' in body,'*' in body,h[1] if h else 'none',h[2] if h else '',h[6] if h else '',round(100*int(h[3])/int(h[5]),2) if h else '',h[7] if h else '',h[8] if h else '',h[9] if h else ''])
for label in ('lt50aa','gt5000aa','contains_X'):print(label,counts[label+'_total'],counts[label+'_with_hit'])
# HSP subject coverage is an approximation based on alignment columns, not interval-union coverage.
