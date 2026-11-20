#!/usr/bin/env python3
"""Draw the user-reported IQ-TREE phylogram with large fonts and a real scale.
Requirements: Python >=3.9, matplotlib. No BioPython needed.
"""
from pathlib import Path
from dataclasses import dataclass, field
import csv, json, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
NEWICK = '(D_latus:0.0542343942,(((((E_granulosus:0.0103170797,E_multilocularis:0.0166552776)100/100:0.0574024437,(T_multiceps:0.0254810903,T_solium:0.0188929201)100/100:0.0592455815)100/100:0.1447749396,H_microstoma:0.3496462277)100/100:0.6968277631,(S_proliferum:0.0251479246,Spirometra_Aus1:0.0203701000)100/100:0.0472178871)100/100:0.0191517158,(S_cotti:0.0102879817,S_solidus:0.0175745143)100/100:0.0419465743)100/100:0.1021299206,L_intestinalis:0.0536221917);'

@dataclass(eq=False)
class Node:
    name: str = ''
    length: float = 0.0
    children: list = field(default_factory=list)

def parse(s):
    i = 0
    def tok():
        nonlocal i
        j = i
        while i < len(s) and s[i] not in ',():;': i += 1
        return s[j:i].strip()
    def sub():
        nonlocal i
        n = Node()
        if s[i] == '(':
            i += 1
            n.children.append(sub())
            while s[i] == ',':
                i += 1
                n.children.append(sub())
            assert s[i] == ')'
            i += 1
        n.name = tok()
        if i < len(s) and s[i] == ':':
            i += 1
            n.length = float(tok())
        return n
    root = sub()
    assert s[i:] == ';'
    return root

def tips(n):
    return {n.name} if not n.children else set().union(*(tips(c) for c in n.children))

def walk(n):
    yield n
    for c in n.children: yield from walk(c)

original = parse(NEWICK)
ALL = tips(original)
C = {'H_microstoma','E_multilocularis','E_granulosus','T_multiceps','T_solium'}
P = {'S_proliferum','Spirometra_Aus1'}
S = {'S_cotti','S_solidus'}
D = {'D_latus','L_intestinalis'}
assert ALL == C | P | S | D

def key(side):
    a, b = tuple(sorted(side)), tuple(sorted(ALL-side))
    return min(a,b)

# Statistics are associated with unrooted edges, not node names.
GCF = {key(C|P):38.61, key(C):99.2, key(P):87.4,
       key(S):86.06, key(D):85.52,
       key({'E_multilocularis','E_granulosus','T_multiceps','T_solium'}):93.57,
       key({'E_multilocularis','E_granulosus'}):95.17,
       key({'T_multiceps','T_solium'}):93.30}
SCF = {key(C|P):31.89, key(C):96.82, key(P):81.29,
       key(S):92.20, key(D):94.16,
       key({'E_multilocularis','E_granulosus','T_multiceps','T_solium'}):81.05,
       key({'E_multilocularis','E_granulosus'}):91.60,
       key({'T_multiceps','T_solium'}):86.32}
SUPPORT = {key(tips(n)):n.name for n in walk(original) if n.children and n is not original}
assert set(GCF) == set(SUPPORT) == set(SCF)

# Reroot on the edge separating the two sampled orders; preserve distances.
adj = {n:[] for n in walk(original)}
for n in walk(original):
    for child in n.children:
        adj[n].append((child, child.length))
        adj[child].append((n, child.length))
clade = next(n for n in walk(original) if tips(n) == C)
parent, edge = next((n,l) for n,l in adj[clade] if tips(n) != C and n not in clade.children)
def orient(node, previous, length):
    return Node(node.name, length, [orient(n,node,l) for n,l in adj[node] if n is not previous])
root = Node('', 0, [orient(clade,parent,edge/2), orient(parent,clade,edge/2)])
preference = ['E_granulosus','E_multilocularis','T_multiceps','T_solium','H_microstoma',
              'S_cotti','S_solidus','D_latus','L_intestinalis','S_proliferum','Spirometra_Aus1']
rank = {s:i for i,s in enumerate(preference)}
def rotate(n):
    for c in n.children: rotate(c)
    n.children.sort(key=lambda c:min(rank[s] for s in tips(c)))
rotate(root)

def distances(r):
    result={}
    def paths(n, path):
        path = dict(path);path[n] = n.length
        if not n.children: result[n.name]=path
        for c in n.children: paths(c,path)
    paths(r,{})
    return {(a,b):sum(v for k,v in result[a].items() if k not in result[b])+sum(v for k,v in result[b].items() if k not in result[a])
            for a in ALL for b in ALL if a < b}
old, new = distances(original), distances(root)
error=max(abs(old[k]-new[k]) for k in old)
assert error < 1e-10

def newick(n):
    body='('+','.join(newick(c) for c in n.children)+')' if n.children else n.name
    label=SUPPORT.get(key(tips(n)),'') if n.children and n is not root else ''
    return body+label+f':{n.length:.10f}'
(HERE/'IQTREE.original.unrooted.nwk').write_text(NEWICK+'\n')
(HERE/'IQTREE.cyclophyllidea_rooted.nwk').write_text(newick(root)+';\n')

labels={'E_granulosus':'E. granulosus','E_multilocularis':'E. multilocularis',
'T_multiceps':'T. multiceps','T_solium':'T. solium','H_microstoma':'H. microstoma',
'S_cotti':'S. cotti','S_solidus':'S. solidus','D_latus':'D. latus',
'L_intestinalis':'L. intestinalis','S_proliferum':'Sparganum proliferum',
'Spirometra_Aus1':'Spirometra sp. Aus1'}
leaves=[n for n in walk(root) if not n.children]
y={n.name:i for i,n in enumerate(leaves)}
xpos={};ypos={}
def layout(n,x=0):
    xpos[n]=x
    for c in n.children:layout(c,x+c.length)
    ypos[n]=sum(ypos[c] for c in n.children)/len(n.children) if n.children else y[n.name]
layout(root)
conflict=key(C|P)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':22,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(22,13.5))
fig.subplots_adjust(left=.045,right=.985,bottom=.24,top=.87)
for sp,col in [('S_cotti','#D97706'),('S_solidus','#C73F38')]:
    ax.axhspan(y[sp]-.4,y[sp]+.4,color=col,alpha=.065,lw=0)
for n in walk(root):
    if n.children:
        ax.plot([xpos[n]]*2,[min(ypos[c] for c in n.children),max(ypos[c] for c in n.children)],color='#343A40',lw=2.2)
        for c in n.children:
            critical=bool(c.children and key(tips(c))==conflict)
            ax.plot([xpos[n],xpos[c]],[ypos[c]]*2,color='#D97706' if critical else '#343A40',lw=4 if critical else 2.2)
            # Annotate each unrooted internal edge once. The root edge is split.
            if c.children and (n is not root or tips(c)==C):
                ax.text(xpos[c]-.004,ypos[c]-.20,SUPPORT[key(tips(c))],ha='right',va='bottom',fontsize=21,
                    color='#AD5D00' if critical else '#343A40',bbox=dict(facecolor='white',edgecolor='none',pad=1.2))
    else:
        col='#008979' if n.name in C else '#86507D'
        if n.name=='S_cotti':col='#D97706'
        if n.name=='S_solidus':col='#C73F38'
        ax.text(xpos[n]+.013,ypos[n],labels[n.name],va='center',fontsize=29,fontstyle='italic',
                fontweight='bold' if n.name in S else 'normal',color=col)
maxx=max(xpos.values())
# Scale is positive and exactly proportional to the drawn branch lengths.
scale=.05
sx=.025;sy=len(leaves)-.05
ax.plot([sx,sx+scale],[sy,sy],color='#252525',lw=3)
ax.plot([sx,sx],[sy-.10,sy+.10],color='#252525',lw=2)
ax.plot([sx+scale,sx+scale],[sy-.10,sy+.10],color='#252525',lw=2)
ax.text(sx+scale/2,sy+.22,'0.05',ha='center',va='top',fontsize=25)
ax.text(sx,sy+.70,'Amino-acid substitutions per site',ha='left',va='top',fontsize=22)
ax.set_xlim(-.012,maxx+.43)
ax.set_ylim(len(leaves)+1.1,-.7)
ax.axis('off')
fig.text(.045,.948,'Single-copy protein phylogeny of eleven cestodes',fontsize=33,fontweight='bold')
fig.text(.045,.905,'373 genes  |  190,437 aligned amino-acid sites  |  Partitioned maximum likelihood',fontsize=23,color='#555555')
fig.legend(handles=[Patch(color='#008979',label='Cyclophyllidea'),Patch(color='#86507D',label='Diphyllobothriidea')],
    loc='lower left',bbox_to_anchor=(.038,.142),ncol=2,frameon=False,fontsize=23)
fig.text(.045,.125,'Branch labels: SH-aLRT (%) / ultrafast bootstrap (%). Root placed between the sampled orders.',fontsize=20,color='#555555')
fig.text(.045,.084,'Orange internal edge: conflicting deep relationship; gCF = 38.61%, sCF = 31.89%.',fontsize=22,color='#AD5D00')
fig.text(.045,.047,'S. cotti + S. solidus: gCF = 86.06%, sCF = 92.20%.',fontsize=22,color='#555555')
for ext in ['png','pdf','svg']:
    fig.savefig(HERE/f'cestode_IQTREE.large_font.{ext}',dpi=250,facecolor='white')
fig.savefig(HERE/'preview.png',dpi=95,facecolor='white')
plt.close(fig)
with (HERE/'branch_support.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['split_side','SH_aLRT_UFBoot','gCF','sCF'])
    for k in GCF:w.writerow([','.join(k),SUPPORT[k],GCF[k],SCF[k]])
(HERE/'validation.json').write_text(json.dumps({'tip_count':len(leaves),'max_pairwise_distance_change':error,'scale_substitutions_per_site':scale,'branch_lengths_proportional':True,'source':'Exact IQ-TREE Newick and concordance values pasted by the user; no server access'},indent=2)+'\n')
(HERE/'README.txt').write_text('Large-font IQ-TREE phylogram. Branch lengths are proportional to amino-acid substitutions per site. Scale = 0.05. Source: user-pasted results. Rerooted between Cyclophyllidea and Diphyllobothriidea; all pairwise tip distances verified unchanged. Supports are associated with unrooted edges. Rooting does not independently establish the root. No gene-duplication counts were moved onto this topology. Orange edge flags concordance conflict; not a low bootstrap branch. To reproduce locally or on server with matplotlib installed: python3 plot_iqtree.py\n')
print(HERE)
print('Validation:',error)
