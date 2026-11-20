#!/usr/bin/env python3
"""Discover both-lineage summaries; preserve source paths and flag duplicate runs.
No automatic quality labels or silent selection among conflicting historical runs.
"""
from pathlib import Path
from collections import Counter
import argparse,csv,re
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
p.add_argument('--species',nargs='+',default=['S_solidus','L_intestinalis','Spirometra_Aus1','D_latus','S_proliferum','E_multilocularis','H_microstoma','E_granulosus','T_solium','T_multiceps'])
a=p.parse_args()
pattern=re.compile(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)')
rows=[]
for path in a.root.rglob('short_summary*.txt'):
 text=path.read_text();match=pattern.search(text)
 if not match:continue
 taxa=[sp for sp in a.species if re.search(r'(?<![A-Za-z0-9_])'+re.escape(sp)+r'(?![A-Za-z0-9_])',str(path))]
 lineages=[lin for lin in ('metazoa_odb10','metazoa_odb12') if lin in str(path) or lin in text]
 if len(taxa)!=1 or len(lineages)!=1:continue
 lin=lineages[0];expected=954 if lin.endswith('10') else 672
 if int(match.group(6))!=expected:raise ValueError(f'Unexpected dataset size: {path}')
 rows.append([taxa[0],lin,*match.groups(),str(path)])
if not rows:raise ValueError('No summaries found; verify --root')
counts=Counter((r[0],r[1]) for r in rows)
a.out.parent.mkdir(parents=True,exist_ok=True)
with a.out.open('x') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['species','lineage','C','S','D','F','M','n','source','duplicate_species_lineage'])
 for r in sorted(rows):w.writerow(r+[counts[(r[0],r[1])]>1])
print('Summary files:',len(rows))
for key,n in sorted(counts.items()):
 if n>1:print('REVIEW duplicate:',*key,n)
for sp in a.species:
 for lin in ('metazoa_odb10','metazoa_odb12'):
  if not counts[(sp,lin)]:print('MISSING:',sp,lin)
print('Output:',a.out)
