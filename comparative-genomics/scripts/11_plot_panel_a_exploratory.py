#!/usr/bin/env python3
"""Reproduce exploratory Panel A from the exported OrthoFinder result subset.
Python >=3.9 and matplotlib required. No topology edits, statistical gain/loss tests,
or automatic TE filtering are performed. Run: python plot_panel_a.py --source source --out figure
"""
from pathlib import Path
from dataclasses import dataclass,field
from collections import Counter
import argparse,csv,hashlib,json,math,re,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Patch
from matplotlib.ticker import FuncFormatter

@dataclass(eq=False)
class Node:
    name:str=''
    length:float=0.0
    children:list=field(default_factory=list)

def read_newick(path):
    s=path.read_text().strip();i=0
    def token():
        nonlocal i
        start=i
        while i<len(s) and s[i] not in ',():;':i+=1
        return s[start:i].strip()
    def parse():
        nonlocal i
        n=Node()
        if s[i]=='(':
            i+=1;n.children=[parse()]
            while s[i]==',':i+=1;n.children.append(parse())
            if s[i]!=')':raise ValueError('Invalid Newick')
            i+=1
        n.name=token()
        if i<len(s) and s[i]==':':i+=1;n.length=float(token())
        return n
    root=parse()
    if s[i:]!=';':raise ValueError('Unexpected Newick suffix')
    return root

def leaves(n):
    return [n.name] if not n.children else [x for c in n.children for x in leaves(c)]
def walk(n):
    yield n
    for child in n.children:yield from walk(child)
def clade_key(n):return frozenset(leaves(n))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    here=Path(__file__).resolve().parent
    p.add_argument('--source',type=Path,default=here/'source')
    p.add_argument('--out',type=Path,default=here/'figure')
    p.add_argument('--minimal',action='store_true',
                   help='Draw only the four data panels; omit overall title, subtitles, legend, notes and caption file')
    p.add_argument('--font-scale',type=float,default=1.0,
                   help='Multiply all plotted font sizes (recommended clean figure: 1.25)')
    a=p.parse_args();src=a.source;out=a.out;out.mkdir(parents=True,exist_ok=True)
    if a.font_scale <= 0:
        raise ValueError('--font-scale must be greater than zero')
    fs=lambda size:size*a.font_scale
    paths={
      'tree':src/'Species_Tree/SpeciesTree_rooted_node_labels.txt',
      'support_tree':src/'Species_Tree/SpeciesTree_rooted.txt',
      'dup_tree':src/'Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt',
      'duplicates':src/'Comparative_Genomics_Statistics/Duplications_per_Species_Tree_Node.tsv',
      'events':src/'Gene_Duplication_Events/Duplications.tsv',
      'per_species':src/'Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv',
      'overall':src/'Comparative_Genomics_Statistics/Statistics_Overall.tsv',
      'gene_counts':src/'Orthogroups/Orthogroups.GeneCount.tsv',
    }
    tree=read_newick(paths['tree']);support_tree=read_newick(paths['support_tree'])
    tree_species=set(leaves(tree));assert len(tree_species)==11
    assert tree_species==set(leaves(support_tree))
    supports={clade_key(n):n.name for n in walk(support_tree) if n.children}
    assert {clade_key(n) for n in walk(tree)}=={clade_key(n) for n in walk(support_tree)},'Tree topologies disagree'
    with paths['per_species'].open() as f:
        rows=list(csv.reader(f,delimiter='\t'))
    species=rows[0][1:];assert set(species)==tree_species
    metrics={}
    for row in rows[1:]:
        if not row or not row[0].strip():break
        if len(row)!=len(species)+1:raise ValueError('Malformed first per-species statistics block')
        metrics[row[0]]={sp:float(v) for sp,v in zip(species,row[1:])}
    with paths['duplicates'].open() as f:
        dup_rows=list(csv.DictReader(f,delimiter='\t'))
    dups={r['Species Tree Node']:int(r['Duplications (50% support)']) for r in dup_rows}
    dups_all={r['Species Tree Node']:int(r['Duplications (all)']) for r in dup_rows}
    assert set(dups)=={n.name for n in walk(tree)}
    event_counts=Counter()
    with paths['events'].open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        for row in reader:
            if float(row['Support'])>=0.5:event_counts[row['Species Tree Node']]+=1
    assert all(event_counts[node]==n for node,n in dups.items()),'Duplication event/table mismatch'
    totals=Counter();unique_og=Counter();unique_genes=Counter();og_n=0
    with paths['gene_counts'].open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        assert set(reader.fieldnames)-{'Orthogroup','Total'}==tree_species
        for row in reader:
            counts={sp:int(row[sp]) for sp in species}
            assert sum(counts.values())==int(row['Total'])
            totals.update(counts);og_n+=1
            present=[sp for sp,n in counts.items() if n]
            if len(present)==1:
                sp=present[0];unique_og[sp]+=1;unique_genes[sp]+=counts[sp]
    required=['Number of genes','Number of genes in orthogroups','Percentage of genes in orthogroups','Number of species-specific orthogroups','Number of genes in species-specific orthogroups','Number of unassigned genes']
    assert all(k in metrics for k in required)
    gene_count_totals = totals.copy()
    count_deltas = {sp: gene_count_totals[sp]-int(metrics['Number of genes in orthogroups'][sp]) for sp in species}
    totals = Counter({sp:int(metrics['Number of genes in orthogroups'][sp]) for sp in species})
    for sp in species:
        assert unique_og[sp]==int(metrics['Number of species-specific orthogroups'][sp])
        assert unique_genes[sp]==int(metrics['Number of genes in species-specific orthogroups'][sp])
        assert metrics['Number of genes'][sp]==totals[sp]+metrics['Number of unassigned genes'][sp]
        actual=100*totals[sp]/metrics['Number of genes'][sp]
        assert abs(actual-metrics['Percentage of genes in orthogroups'][sp])<=0.051
    overall={}
    for line in paths['overall'].read_text().splitlines():
        if not line.strip():break
        f=line.split('\t');overall[f[0]]=f[1]
    statistics_og_n=int(overall['Number of orthogroups'])
    assert statistics_og_n==16669
    assert sum(totals.values())==int(overall['Number of genes in orthogroups'])==137591
    assert sum(metrics['Number of genes'].values())==int(overall['Number of genes'])==161409
    preference=['T_multiceps','T_solium','E_granulosus','E_multilocularis','H_microstoma','S_solidus','S_cotti','Spirometra_Aus1','S_proliferum','D_latus','L_intestinalis']
    rank={sp:i for i,sp in enumerate(preference)}
    def rotate(n):
        for c in n.children:rotate(c)
        n.children.sort(key=lambda c:min(rank[x] for x in leaves(c)))
    before={clade_key(n) for n in walk(tree)};rotate(tree)
    assert before=={clade_key(n) for n in walk(tree)},'Rotation changed topology'
    order=leaves(tree);y={sp:i for i,sp in enumerate(order)}
    def height(n):return 0 if not n.children else 1+max(height(c) for c in n.children)
    H=height(tree);xx={n:H-height(n) for n in walk(tree)}
    yy={}
    def set_y(n):
        if not n.children:yy[n]=y[n.name]
        else:
            for c in n.children:set_y(c)
            yy[n]=sum(yy[c] for c in n.children)/len(n.children)
    set_y(tree)
    cyclo={'T_multiceps','T_solium','E_granulosus','E_multilocularis','H_microstoma'}
    diphy=tree_species-cyclo
    diphy_monophyletic=any(clade_key(n)==diphy for n in walk(tree))
    names={'T_multiceps':'T. multiceps','T_solium':'T. solium','E_granulosus':'E. granulosus','E_multilocularis':'E. multilocularis','H_microstoma':'H. microstoma','S_solidus':'S. solidus','S_cotti':'S. cotti','Spirometra_Aus1':'Spirometra sp. Aus1','S_proliferum':'Sparganum proliferum','D_latus':'D. latus','L_intestinalis':'L. intestinalis'}
    colors={sp:('#258777' if sp in cyclo else '#865277') for sp in species}
    colors.update(S_solidus='#b65147',S_cotti='#d4853c')
    with (out/'panel_A_species_data.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['species','display_name','order','total_input_genes','terminal_duplications_support_ge_0.5','genes_in_orthogroups','percent_in_orthogroups','species_specific_orthogroups','genes_in_species_specific_orthogroups','unassigned_genes'])
        for sp in order:w.writerow([sp,names[sp],'Cyclophyllidea' if sp in cyclo else 'Diphyllobothriidea',int(metrics['Number of genes'][sp]),dups[sp],totals[sp],metrics['Percentage of genes in orthogroups'][sp],unique_og[sp],unique_genes[sp],int(metrics['Number of unassigned genes'][sp])])
    with (out/'panel_A_node_data.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['node','descendant_species','duplications_all','duplications_support_ge_0.5','exported_tree_support'])
        for n in walk(tree):
            if n.children:w.writerow([n.name,','.join(leaves(n)),dups_all[n.name],dups[n.name],supports[clade_key(n)]])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':fs(10),'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.65})
    fig=plt.figure(figsize=(18.8,7.2 if a.minimal else 7.7),facecolor='white')
    bottom=.105 if a.minimal else .205
    top=.855 if a.minimal else .78
    gs=fig.add_gridspec(1,5,width_ratios=[5.25,.15,2.6,2.85,2.75],left=.025,right=.985,bottom=bottom,top=top,wspace=.20)
    tree_ax=fig.add_subplot(gs[0]);strip=fig.add_subplot(gs[1]);axes=[fig.add_subplot(gs[i]) for i in (2,3,4)]
    all_axes=[tree_ax,strip,*axes]
    for ax in all_axes:
        ax.set_ylim(len(order)-.4,-.65)
        for sp in ['S_solidus','S_cotti']:
            ax.axhspan(y[sp]-.43,y[sp]+.43,color=colors[sp],alpha=.075,zorder=0,lw=0)
    # Draw a topology-only cladogram. Branch lengths are intentionally not represented.
    for n in walk(tree):
        if n.children:
            tree_ax.plot([xx[n],xx[n]],[min(yy[c] for c in n.children),max(yy[c] for c in n.children)],color='#42464b',lw=1.15,zorder=2)
            for c in n.children:tree_ax.plot([xx[n],xx[c]],[yy[c],yy[c]],color='#42464b',lw=1.15,zorder=2)
            support=supports[clade_key(n)]
            label=f'{dups[n.name]:,}' + (f' ({float(support):g})' if support else '')
            tree_ax.text(xx[n]-.05,yy[n]-.19,label,ha='right',va='bottom',fontsize=fs(8.8),color='#4c5056',bbox=dict(facecolor='white',edgecolor='none',pad=.25),zorder=4)
        else:
            tree_ax.text(H+.18,yy[n],names[n.name],va='center',fontsize=fs(12.2),fontstyle='italic',fontweight='bold' if n.name in {'S_solidus','S_cotti'} else 'normal',color=colors[n.name])
    tree_ax.set_xlim(-1.05,H+4.05);tree_ax.axis('off')
    tree_ax.text(0,1.08 if a.minimal else 1.13,'Species phylogeny',transform=tree_ax.transAxes,fontsize=fs(14.2),fontweight='bold')
    if not a.minimal:
        tree_ax.text(0,1.065,'Internal duplications (exported tree support)',transform=tree_ax.transAxes,fontsize=fs(9.2),color='#62676c')
    strip.set_xlim(0,1);strip.axis('off')
    for sp in order:strip.add_patch(Rectangle((0,y[sp]-.43),1,.86,color=colors[sp] if sp not in {'S_solidus','S_cotti'} else '#865277',lw=0))
    headings=['Terminal gene\nduplications','Proteins in\northogroups','Species-specific\northogroups']
    subtitles=['Events; support ≥0.5','Bar = assigned %; label = gene count','Groups represented in one species only']
    values=[[dups[sp] for sp in order],[metrics['Percentage of genes in orthogroups'][sp] for sp in order],[unique_og[sp] for sp in order]]
    xmax=[14200,128,1830];ticks=[[0,5000,10000],[0,25,50,75,100],[0,500,1000,1500]]
    for j,ax in enumerate(axes):
        vals=values[j]
        ax.barh(range(len(order)),vals,height=.68,color=[colors[sp] for sp in order],edgecolor='none',zorder=3)
        ax.set_xlim(0,xmax[j]);ax.set_xticks(ticks[j]);ax.set_yticks([]);ax.xaxis.grid(True,color='#e8ebee',lw=.7,zorder=1);ax.set_axisbelow(True)
        for side in ['top','right','left']:ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color('#c4c9cd');ax.tick_params(axis='x',length=3,color='#a0a7ad',labelsize=fs(9.5))
        if j!=1:ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{int(x):,}'))
        ax.set_xlabel('Assigned proteins (%)' if j==1 else ('Events' if j==0 else 'Orthogroups'),labelpad=8,fontsize=fs(10.2))
        ax.text(0,1.08 if a.minimal else 1.13,headings[j],transform=ax.transAxes,fontsize=fs(13.2),fontweight='bold',linespacing=1.05)
        if not a.minimal:
            ax.text(0,1.065,subtitles[j],transform=ax.transAxes,fontsize=fs(8.6),color='#62676c')
        for i,(sp,value) in enumerate(zip(order,vals)):
            label=f'{totals[sp]:,}' if j==1 else f'{int(value):,}'
            ax.text(value+xmax[j]*.018,i,label,ha='left',va='center',fontsize=fs(10.2),fontweight='bold' if sp in {'S_solidus','S_cotti'} else 'normal',color='#30343a',zorder=5)
    if not a.minimal:
        fig.text(.027,.945,'a',fontsize=fs(22),fontweight='bold')
        fig.text(.065,.947,'Orthology and inferred gene duplications across eleven cestodes',fontsize=fs(17),fontweight='bold',color='#26343c')
        fig.text(.065,.905,'EXPLORATORY RUN  •  cotti filtered7 candidate  •  TE-associated proteins not uniformly filtered',fontsize=fs(10),color='#72777d')
        legend=[Patch(color='#258777',label='Cyclophyllidea'),Patch(color='#865277',label='Diphyllobothriidea'),Patch(color='#b65147',label='S. solidus'),Patch(color='#d4853c',label='S. cotti (candidate)')]
        fig.legend(handles=legend,loc='lower left',bbox_to_anchor=(.06,.105),ncol=4,frameon=False,fontsize=fs(10),handlelength=1.1,columnspacing=2.2)
        fig.text(.065,.075,'Topology is shown as a cladogram (branch lengths not to scale); node supports are reproduced as exported, not called bootstraps.',fontsize=fs(9),color='#656a70')
        fig.text(.065,.047,'Order colours denote taxonomy, not monophyly. The inferred tree does not recover Diphyllobothriidea as monophyletic; topology review is required.',fontsize=fs(9),color='#656a70')
    output_stem='panel_A_clean' if a.minimal else 'panel_A_exploratory'
    for ext in ['pdf','svg','png']:
        fig.savefig(out/f'{output_stem}.{ext}',dpi=300,facecolor='white',bbox_inches='tight' if a.minimal else None)
    fig.savefig(out/f'{output_stem}_preview.png',dpi=130,facecolor='white',bbox_inches='tight' if a.minimal else None);plt.close(fig)
    caption='''Panel A. Exploratory orthology and gene-duplication summary for eleven cestodes. The left panel reproduces the rooted OrthoFinder topology as a cladogram; branch lengths are not represented. Internal labels give inferred duplication events with duplication support at least 0.5, followed in parentheses by tree support values exactly as exported by OrthoFinder (not interpreted as bootstrap percentages). The root has no exported support value. Terminal duplication counts use the same event-support threshold. The middle bars show the percentage of input representative proteins assigned to orthogroups; labels show assigned gene counts. The right bars show orthogroups represented in only one of the eleven input species; unassigned singleton genes are not included in this count. Colours denote taxonomic order, with the two focal species highlighted. The cotti input is the filtered7 working candidate. TE-associated proteins have not been uniformly removed, and annotation completeness differs among taxa. Accordingly, species-specific groups and duplication counts are descriptive outputs rather than validated lineage-specific innovations or family expansions. The supplied rooted tree places D. latus and L. intestinalis outside the other nine species and does not recover Diphyllobothriidea as monophyletic; this topology is retained transparently and requires independent review before final phylogenetic interpretation.\n'''
    caption += '\nData-source audit: assigned gene counts and percentages in the figure use Statistics_PerSpecies.tsv. Orthogroups.GeneCount.tsv contains 16,583 rows and 137,611 summed memberships, whereas Statistics_Overall.tsv reports 16,669 orthogroups and 137,591 assigned genes. The reason is unresolved in this exported subset; the tables must not be treated as interchangeable. Species-specific group and gene counts agree exactly between the two sources. No assumption about the reason for this discrepancy was made.\n'
    if not a.minimal:
        (out/'figure_caption.txt').write_text(caption)
    validation=dict(species=11,input_genes=161409,genes_in_orthogroups=137591,statistics_orthogroups=statistics_og_n,gene_count_table_orthogroups=og_n,gene_count_table_memberships=sum(gene_count_totals.values()),gene_count_minus_statistics_by_species=count_deltas,single_copy_orthogroups=int(overall['Number of single-copy orthogroups']),duplication_event_threshold=0.5,duplication_counts_match_event_table=True,gene_counts_match_statistics=not any(count_deltas.values()),species_specific_counts_match_gene_count_table=True,topology_preserved=True,diphyllobothriidea_monophyletic=diphy_monophyletic,tip_order=order,matplotlib_version=matplotlib.__version__,python_version=sys.version)
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    with (out/'input.SHA256SUMS').open('w') as f:
        for label,path in paths.items():f.write(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n')
    print(json.dumps(validation,indent=2))
    print('Figure outputs:',out)

if __name__=='__main__':main()
