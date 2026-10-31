#!/usr/bin/env python3
"""Build a review candidate from the original S. cotti MITOS2 GFF, not a submission file.
Requires Biopython in the remote mitos2 environment. Never edits the source FASTA/GFF.
Circular origin-spanning RNA features retain GFF3 virtual coordinates.
"""
import argparse
import csv
import hashlib
import itertools
import re
import shutil
from collections import Counter
from pathlib import Path
from urllib.parse import quote, unquote
from Bio import SeqIO
from Bio.Seq import Seq

PCG = {'nad1','nad2','atp6','nad4','nad4l','cob','cox3','nad5','nad6','cox2','cox1','nad3'}
TRNA = {'trnN','trnD','trnA','trnV','trnM','trnF','trnQ','trnG','trnR','trnL2','trnS2','trnL1','trnY','trnE','trnC','trnT','trnW','trnS1','trnK','trnI','trnP','trnH'}
BASELINE = {'atp6':(2072,2614), 'cox3':(5366,6016), 'nad3':(14857,15228)}
CHANGES = {
    'atp6': (2072,2587,'Published MW602521.1 full-query homology; 171 aa plus complete stop'),
    'cox3': (5428,6070,'Published MW602521.1 full-query homology; GTG start; 214 aa plus terminal T'),
    'trnH': (5365,5427,'Published MW602521.1 full-query homology; ARAGORN candidate 5364-5428; boundary candidate'),
    'rrnL': (12162,13125,'Published MW602521.1 full-query homology; 964 bp; boundary candidate'),
    'nad3': (14883,15228,'Partial published homology plus terminal review; ATG start; 115 aa plus terminal T'),
}

def attrs(text):
    return dict((k,unquote(v)) for k,v in (item.split('=',1) for item in text.split(';') if '=' in item))

def segments(start,end,length):
    return [(start,end)] if end<=length else [(start,length),(1,end-length)]

def build(gff_path, fasta_path, out):
    records=list(SeqIO.parse(fasta_path,'fasta'))
    if len(records)!=1 or records[0].id!='contig_1' or len(records[0])!=15437:
        raise ValueError('Expected the original contig_1, 15437 bp; do not use a rotated FASTA')
    rec=records[0]; length=len(rec); genes={}; exons={}; raw=[]
    with open(gff_path) as handle:
        for number,line in enumerate(handle,1):
            if line.startswith('##FASTA'): break
            if line.startswith('#') or not line.strip(): continue
            p=line.rstrip('\r\n').split('\t')
            if len(p)!=9: raise ValueError(f'Line {number}: expected 9 TAB-separated columns')
            a=attrs(p[8]); name=a.get('Name',a.get('gene_id','')).split('(')[0]
            raw.append((line,p,a,name))
            if p[0]!=rec.id: raise ValueError(f'Unexpected sequence ID: {p[0]}')
            if p[2] in ('gene','ncRNA_gene') and name in PCG|TRNA|{'rrnL','rrnS'}:
                if name in genes: raise ValueError(f'Duplicate gene: {name}')
                genes[name]=(int(p[3]),int(p[4]),p[6],p[1])
            if p[2]=='exon' and name in PCG:
                exons.setdefault(name,[]).append((int(p[3]),int(p[4]),p[7]))
    expected=PCG|(TRNA-{'trnH'})|{'rrnS'}
    if set(genes)!=expected:
        raise ValueError(f'Unexpected original gene set; missing={expected-set(genes)}, extra={set(genes)-expected}')
    for name in PCG:
        start,end,strand,source=genes[name]
        if exons.get(name)!=[(start,end,'0')]:
            raise ValueError(f'{name}: expected one phase-0 MITOS exon matching the gene; manual review required')
    for name,coords in BASELINE.items():
        if genes[name][:2]!=coords: raise ValueError(f'{name}: unexpected original coordinates')
    original=dict(genes)
    for name,(start,end,evidence) in CHANGES.items():
        genes[name]=(start,end,'-','manual_review')
    if set(genes)!=PCG|TRNA|{'rrnL','rrnS'}: raise ValueError('Unexpected revised gene set')
    rows=[]; proteins=[]; audits=[]
    def add(source,kind,start,end,strand,attributes,phase='.'):
        if not (1<=start<=length and start<=end<start+length):
            raise ValueError(f'Invalid circular coordinates: {kind} {start}-{end}')
        value=';'.join(k+'='+quote(str(v),safe=':@|+-_.') for k,v in attributes.items())
        rows.append([rec.id,source,kind,str(start),str(end),'.',strand,phase,value])
    add('manual_review','region',1,length,'+',{'ID':'Scotti_mt_region','Name':rec.id,'Is_circular':'true'})
    for name,(start,end,strand,source) in sorted(genes.items(),key=lambda item:item[1][0]):
        gid=f'Scotti_mt_gene_{name}'; tid=f'Scotti_mt_transcript_{name}'
        note=CHANGES[name][2] if name in CHANGES else 'Retained MITOS2 prediction; not independently revalidated'
        if end>length: note+='; origin-spanning virtual coordinates retained; boundary requires review'
        if name in {'cox1','trnT'}: note+='; cox1-trnT 10 bp overlap remains under review'
        add(source,'gene' if name in PCG else 'ncRNA_gene',start,end,strand,{'ID':gid,'Name':name,'Note':note,'curation_status':'candidate'})
        add(source,'mRNA' if name in PCG else ('rRNA' if name in {'rrnL','rrnS'} else 'tRNA'),start,end,strand,{'ID':tid,'Parent':gid,'Name':name})
        add(source,'exon',start,end,strand,{'ID':f'{tid}_exon','Parent':tid})
        if name not in PCG: continue
        if end>length: raise ValueError('This script does not support origin-spanning CDS')
        dna=rec.seq[start-1:end]
        if strand=='-': dna=dna.reverse_complement()
        dna=str(dna).upper(); full=len(dna)//3*3; remainder=dna[full:]
        protein=str(Seq(dna[:full]).translate(table=9))
        internal='*' in protein.rstrip('*'); complete=protein.endswith('*') and not remainder
        protein=protein[:-1] if protein.endswith('*') else protein
        if dna[:3] in {'ATG','GTG'} and protein: protein='M'+protein[1:]
        if name in {'cox3','nad3'} and (remainder!='T' or internal or len(protein)!={'cox3':214,'nad3':115}[name]):
            raise ValueError(f'{name}: sequence does not match the verified candidate')
        if name=='atp6' and (internal or not complete or len(protein)!=171):
            raise ValueError('atp6: sequence does not match the verified candidate')
        ca={'ID':f'Scotti_mt_CDS_{name}','Parent':tid,'Name':name,'transl_table':'9'}
        if remainder in {'T','TA'}:
            ca['Note']='Incomplete terminal '+remainder+'; completion by polyadenylation inferred from homology; no genomic bases added'
            ca['incomplete_stop_codon']=remainder
        add(source,'CDS',start,end,strand,ca,'0')
        proteins.append((name,protein))
        audits.append([name,start,end,strand,len(dna),dna[:3],len(protein),remainder or '.',str(complete),str(internal)])
    # Validate rebuilt ID/Parent graph before creating output files.
    ids={attrs(r[8])['ID']:r for r in rows}
    if len(ids)!=len(rows): raise ValueError('Duplicate output IDs')
    for r in rows:
        a=attrs(r[8])
        if 'Parent' in a:
            parent=ids.get(a['Parent'])
            if parent is None or not (int(parent[3])<=int(r[3])<=int(r[4])<=int(parent[4])):
                raise ValueError(f'Invalid parent relation: {a}')
    out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(gff_path,out/'MITOS2.original.gff')
    def table(filename,header,data):
        with open(out/filename,'w',newline='') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(header); w.writerows(data)
    with open(out/'Scotti_mitochondrial.curated_candidate.gff3','w') as f:
        f.write(f'##gff-version 3\n##sequence-region {rec.id} 1 {length}\n')
        f.write('# Review candidate only; circular virtual coordinates are intentional. Not a validated submission.\n')
        for row in rows: f.write('\t'.join(row)+'\n')
    with open(out/'Scotti_mitochondrial.candidate.proteins.fa','w') as f:
        for name,protein in proteins: f.write(f'>{name}\n{protein}\n')
    table('translation_audit.tsv',['gene','start','end','strand','nt_length','start_codon','aa_length','trailing_bases','complete_stop','internal_stop'],audits)
    table('changes.tsv',['gene','old_start','old_end','new_start','new_end','strand','evidence'],[
        [name,*(original[name][:2] if name in original else ('.','.')),start,end,'-',evidence]
        for name,(start,end,evidence) in CHANGES.items()])
    overlaps=[]
    for (a,ga),(b,gb) in itertools.combinations(genes.items(),2):
        n=sum(max(0,min(xe,ye)-max(xs,ys)+1) for xs,xe in segments(*ga[:2],length) for ys,ye in segments(*gb[:2],length))
        if n: overlaps.append([a,b,n,'review_not_automatically_trimmed'])
    table('overlap_review.tsv',['gene1','gene2','overlap_bp','status'],overlaps)
    with open(out/'replication_origin_predictions.unreviewed.tsv','w') as f:
        f.write('seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n')
        for line,p,a,name in raw:
            if p[2]=='origin_of_replication' or re.search(r'(^|[;=_-])(OH|OL)([;=_-]|$)',p[8]): f.write('\t'.join(p)+'\n')
    with open(out/'README.txt','w') as f:
        f.write('Candidate GFF only. Original input preserved. FASTA unchanged.\n')
        f.write('36 genes: 12 PCG, 22 tRNA, 2 rRNA. Other retained MITOS predictions are not independently confirmed.\n')
        f.write('trnP 15381-15445 spans the origin (15381-15437 plus 1-8); allowed circular GFF3 virtual coordinates.\n')
        f.write('Some linear-coordinate GFF tools will not accept circular virtual coordinates; adapt before using those tools.\n')
        f.write('Review overlap_review.tsv and translation_audit.tsv before finalizing. No OH/OL functional assignment.\n')
        f.write('FASTA SHA256: '+hashlib.sha256(fasta_path.read_bytes()).hexdigest()+'\n')
        f.write('Original GFF SHA256: '+hashlib.sha256(gff_path.read_bytes()).hexdigest()+'\n')
    print('Created:',out/'Scotti_mitochondrial.curated_candidate.gff3')
    print('Genes: 36; PCG: 12; tRNA: 22; rRNA: 2. ID/Parent and circular-coordinate checks passed.')
    print('Remaining gene overlaps (not automatically errors):')
    for row in overlaps: print(*row,sep='\t')
    warnings=[r[0] for r in audits if r[-1]=='True']
    print('Genes with internal translation stops:',','.join(warnings) or 'none')
    print('Review candidate only; external validator/submission checks have not been run.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    root=Path('/work/cyu/assembly/worm/scotti/mitochondria')
    p.add_argument('--gff',type=Path,default=root/'05_mitos2/result.gff')
    p.add_argument('--fasta',type=Path,default=root/'03_flye_hifi_random700/assembly/assembly.fasta')
    p.add_argument('--out',type=Path,default=root/'06_curated_candidate')
    a=p.parse_args(); build(a.gff,a.fasta,a.out)
