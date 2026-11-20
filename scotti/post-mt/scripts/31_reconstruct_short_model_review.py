#!/usr/bin/env python3
"""Reconstruct short-model evidence/context review. Never removes gene models.
Without --rna, leaves RNA/masking fields unassessed. --rna scans only 44 short-gene loci.
Requires Biopython only for --rna. Run to a new output directory.
"""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,subprocess
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,default=Path('/work/cyu/annotation/scotti-annotation'))
p.add_argument('--out',type=Path,required=True);p.add_argument('--rna',action='store_true');a=p.parse_args()
qc=a.root/'results/annotation_qc';run=a.root/'results/braker3_etp_v1_retry2'
with (qc/'06_protein_support_audit/protein_support.tsv').open() as f:
    support={r['protein_id']:r for r in csv.DictReader(f,delimiter='\t')}
short={pid:r for pid,r in support.items() if int(r['length_aa'])<50}
if len(short)!=46:raise ValueError('Expected the documented 46 short isoforms')
parents={};genes={};cds=defaultdict(list);by_gene=defaultdict(list)
for line in (run/'braker.gff3').open():
    if line.startswith('##FASTA'):break
    if line.startswith('#'):continue
    f=line.rstrip().split('\t')
    if len(f)!=9:continue
    at=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
    if f[2]=='gene':genes[at['ID']]=(f[0],int(f[3])-1,int(f[4]),f[6])
    elif f[2]=='mRNA':parents[at['ID']]=at['Parent'];by_gene[at['Parent']].append(at['ID'])
    elif f[2]=='CDS':
        for tx in at.get('Parent','').split(','):cds[tx].append((f[0],int(f[3])-1,int(f[4])))
egg={};header=None
for line in (a.root/'results/functional/02_eggnog/scotti_braker.emapper.annotations').open():
    if not line.strip():continue
    if line.startswith('#'):
        f=line.lstrip('#').rstrip('\n').split('\t')
        if f[0] in ('query','query_name'):header=f
        continue
    if header is None:raise ValueError('eggNOG header missing')
    f=line.rstrip('\n').split('\t')
    if f[0] in short:egg[f[0]]=dict(zip(header,f))
signatures=defaultdict(set);ipr=defaultdict(set)
for line in (a.root/'results/functional/03_interproscan/scotti_braker.tsv').open():
    if line.startswith('#'):continue
    f=line.rstrip('\n').split('\t')
    if len(f)<11 or f[0] not in short:continue
    signatures[f[0]].add(f'{f[3]}:{f[4]}:{f[5]}')
    if len(f)>=13 and f[11].startswith('IPR'):ipr[f[0]].add(f'{f[11]}:{f[12]}')
genome=None
if a.rna:
    from Bio import SeqIO
    genome=SeqIO.index('/work/cyu/assembly/worm/scotti/repeats/02_repeatmasker/Scotti_nuclear.mt_filtered.ragtag.provisional.fa.masked','fasta')
a.out.mkdir(parents=True,exist_ok=False);counts=Counter();candidate_genes=set()
fields=['protein_id','gene_id','length_aa','longest_isoform_aa','gene_context','swissprot_hit','eggnog_row_present','member_signatures','interpro_entries','CDS_softmasked_pct','unique_primary_records_in_gene_span','review_category','final_decision']
with (a.out/'short_model_review.tsv').open('x') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
    for pid,r in sorted(short.items()):
        gene=parents[pid];maximum=max(int(support[tx]['length_aa']) for tx in by_gene[gene]);length=int(r['length_aa'])
        functional=r['swissprot_hit'] not in ('','none') or pid in egg or bool(signatures[pid]) or bool(ipr[pid])
        context='all_isoforms_below_50aa' if maximum<50 else 'has_isoform_at_least_50aa'
        if maximum>=50:category='short_isoform_in_longer_gene'
        elif length<=10 and not functional:category='priority_exclusion_candidate';candidate_genes.add(gene)
        elif functional:category='retain_pending_review'
        else:category='manual_review'
        mask='not_assessed';rna='not_assessed'
        if a.rna and maximum<50:
            seq=''.join(str(genome[c].seq[s:e]) for c,s,e in cds[pid])
            if not seq:raise ValueError('No CDS sequence')
            mask=round(100*sum(c in 'acgt' for c in seq)/len(seq),2)
            chrom,start,end,strand=genes[gene]
            result=subprocess.run(['samtools','view','-c','-F','2308','-q','255',str(a.root/'results/star/02_rna_alignment/Scotti_nuclear.rnaseq.bam'),f'{chrom}:{start+1}-{end}'],check=True,capture_output=True,text=True)
            rna=int(result.stdout.strip())
        counts[category]+=1
        w.writerow(dict(protein_id=pid,gene_id=gene,length_aa=length,longest_isoform_aa=maximum,gene_context=context,swissprot_hit=r['swissprot_hit'],eggnog_row_present=pid in egg,member_signatures=';'.join(sorted(signatures[pid])),interpro_entries=';'.join(sorted(ipr[pid])),CDS_softmasked_pct=mask,unique_primary_records_in_gene_span=rna,review_category=category,final_decision='not_finalized'))
if genome is not None:genome.close()
(a.out/'priority_candidate_gene_ids.txt').write_text(''.join(x+'\n' for x in sorted(candidate_genes)))
for category,n in sorted(counts.items()):print(category,n)
print('No annotation sequences were removed. Final seven-model exclusion is documented separately.')
