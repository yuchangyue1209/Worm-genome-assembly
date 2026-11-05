#!/usr/bin/env python3
"""Compare odb10 genome/all/longest statuses and reproduce Miniprot stop metadata.
Requires a new --out directory. Does not modify predictions or rerun BUSCO.
"""
from pathlib import Path
from collections import defaultdict,Counter
from Bio import SeqIO
import argparse,csv
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
root=Path('/work/cyu/annotation/scotti-annotation')
base=Path('/work/cyu/assembly/worm/scotti/qc/primary_busco/primary_metazoa_odb10')
run=base/'run_metazoa_odb10'
paths={'genome':run/'full_table.tsv','all_proteins':root/'results/annotation_qc/01_raw_braker_busco/raw_braker_metazoa/run_metazoa_odb10/full_table.tsv','longest':root/'results/annotation_qc/02_longest_isoform_busco/longest_metazoa/run_metazoa_odb10/full_table.tsv'}
states={}
for label,path in paths.items():
    d={}
    for line in path.open():
        if not line.strip() or line.startswith('#'): continue
        f=line.rstrip().split('	')
        if f[0] in d and d[f[0]]!=f[1]: raise ValueError('Inconsistent BUSCO statuses')
        d[f[0]]=f[1]
    states[label]=d
ids=set(states['genome'])
assert len(ids)==954 and all(set(x)==ids for x in states.values())
selected={};wanted=set()
for folder in ('single_copy_busco_sequences','multi_copy_busco_sequences'):
    for path in (run/'busco_sequences'/folder).glob('*.faa'):
        names={r.id for r in SeqIO.parse(path,'fasta')}
        if not names: raise ValueError(path)
        selected[path.stem]=names;wanted.update(names)
metadata={}
for line in (base/'logs/miniprot_align_metazoa_odb10_out.log').open():
    if not line.startswith('##PAF	'):continue
    f=line.rstrip().split('	')[1:]
    if len(f)<18 or f[5]=='*':continue
    key=f'{f[0]}|{f[5]}:{f[7]}-{f[8]}|{f[4]}'
    if key not in wanted:continue
    tag=f[17].split(':',2)
    if len(tag)!=3 or tag[1]!='i':raise ValueError('Unexpected BUSCO 6.0 stop-count field')
    n=int(tag[2])
    if key in metadata and metadata[key]!=n:raise ValueError('Conflicting metadata')
    metadata[key]=n
if wanted-metadata.keys():raise ValueError('Selected predictions missing metadata')
audit={bid:[metadata[key] for key in sorted(names)] for bid,names in selected.items()}
flagged=sum(all(n>0 for n in x) for x in audit.values())
if len(audit)!=697 or flagged!=216:raise ValueError(f'Unexpected summary: {len(audit)}, {flagged}')
a.out.mkdir(parents=True,exist_ok=False)
complete={'Complete','Duplicated'};counts=Counter()
with (a.out/'status_comparison.tsv').open('x') as f,(a.out/'candidate_miniprot_stop_audit.tsv').open('x') as g:
    w=csv.writer(f,delimiter='	');v=csv.writer(g,delimiter='	')
    w.writerow(['busco_id',*paths]);v.writerow(['busco_id','protein_status','stop_counts','all_selected_copies_have_stops'])
    for bid in sorted(ids):
        row=[bid]+[states[k][bid] for k in paths];w.writerow(row)
        counts[(row[1],row[2])]+=1
        if row[1] in complete and row[2] not in complete:
            ns=audit[bid];v.writerow([bid,row[2],','.join(map(str,ns)),all(n>0 for n in ns)])
for pair,n in sorted(counts.items()):print(*pair,n,sep='	')
print('Validated complete/stop groups:',len(audit),flagged)
