#!/usr/bin/env python3
"""Reconstructed QC and longest-protein selection. Biopython required.
Run with a NEW output directory. Does not replace original annotations.
"""
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import argparse,csv
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--run',type=Path,default=Path('/work/cyu/annotation/scotti-annotation/results/braker3_etp_v1_retry2'))
p.add_argument('--out',type=Path,required=True)
a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=False)
genes=set(); parents={}
for line in (a.run/'braker.gff3').open():
    if line.startswith('##FASTA'): break
    if line.startswith('#'): continue
    f=line.rstrip().split('	')
    if len(f)!=9: continue
    attrs=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
    if f[2]=='gene': genes.add(attrs['ID'])
    elif f[2]=='mRNA':
        tx=attrs['ID']; parent=attrs['Parent']
        if tx in parents or ',' in parent: raise ValueError('Unexpected transcript mapping')
        parents[tx]=parent
proteins={}
for r in SeqIO.parse(a.run/'braker.aa','fasta'):
    if r.id in proteins: raise ValueError('Duplicate protein ID')
    proteins[r.id]=r
if set(proteins)!=set(parents): raise ValueError('Protein/mRNA ID mismatch')
by_gene=defaultdict(list)
for tx,gene in parents.items(): by_gene[gene].append(proteins[tx])
if set(by_gene)!=genes: raise ValueError('Gene/parent mismatch')
def length(r): return len(str(r.seq).rstrip('*'))
selected=[]
with (a.out/'longest_isoform_selection.tsv').open('x') as f:
    w=csv.writer(f,delimiter='	');w.writerow(['gene_id','transcript_id','protein_length_aa','isoform_count'])
    for gene,records in sorted(by_gene.items()):
        best=min(records,key=lambda r:(-length(r),r.id)); selected.append(best)
        w.writerow([gene,best.id,length(best),len(records)])
SeqIO.write(selected,a.out/'Scotti_braker.longest.proteins.fa','fasta')
with (a.out/'protein_qc.tsv').open('x') as f:
    w=csv.writer(f,delimiter='	');w.writerow(['set','proteins','internal_stop','contains_X','lt50aa','lt100aa','median_aa','max_aa'])
    for label,records in [('all',list(proteins.values())),('longest',selected)]:
        bodies=[str(r.seq).upper().removesuffix('*') for r in records] if hasattr(str,'removesuffix') else [str(r.seq).upper()[:-1] if str(r.seq).endswith('*') else str(r.seq).upper() for r in records]
        lengths=sorted(map(len,bodies))
        row=[label,len(bodies),sum('*' in x for x in bodies),sum('X' in x for x in bodies),sum(n<50 for n in lengths),sum(n<100 for n in lengths),lengths[len(lengths)//2],max(lengths)]
        w.writerow(row);print(*row,sep='	')
