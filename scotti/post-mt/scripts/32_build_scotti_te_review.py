#!/usr/bin/env python3
"""Build a reversible TE-protein review table; never remove gene models.
Standard library only. All regexes are screening rules, not TE confirmations.
"""
from pathlib import Path
from collections import defaultdict, Counter
from bisect import bisect_right
from urllib.parse import unquote
import argparse, csv, hashlib, json, re, shutil

PATTERNS = {
    'explicit_TE': r'transposase|retrotranspos(?:on|able)|retroelement|transposable.element|retrovirus|retroviral|\bgag[ -]pol\b|\bpol polyprotein\b',
    'RT_or_integrase': r'reverse.transcriptase|RNA.directed.DNA.polymerase|RNA.dependent.DNA.polymerase|\bintegrase\b|\bRVT_[123]\b|\brve\b|\bRVP\b',
    'generic_associated': r'\bribonuclease H\b|\bRNase[ _-]?H\b|\bDDE\b|\bD[D/E]D\b|\bretroviral aspartyl protease\b|\bAP_endonuc_2\b',
    'host_context': r'telomerase|\bTERT\b|domesticat|retrotransposon.derived|transposase.derived|\bRAG[12]\b|recombination.activating|\bPEG10\b|\bZBED\w*\b|zinc.finger.BED',
}
COMPILED = {k:re.compile(v,re.I) for k,v in PATTERNS.items()}
TE_CLASSES = {'LINE','SINE','LTR','DNA','RC','Retroposon','PLE','Penelope','Helitron','Maverick','Crypton','TIR'}

def merge(intervals):
    result=[]
    for start,end in sorted(intervals):
        if end <= start: raise ValueError('Nonpositive interval')
        if result and start <= result[-1][1]:
            result[-1]=(result[-1][0],max(result[-1][1],end))
        else: result.append((start,end))
    return result

class Coverage:
    def __init__(self, intervals):
        self.intervals=merge(intervals)
        self.starts=[s for s,e in self.intervals]
        self.ends=[e for s,e in self.intervals]
        self.prefix=[0]
        for s,e in self.intervals:self.prefix.append(self.prefix[-1]+e-s)
    def before(self,x):
        i=bisect_right(self.ends,x)
        answer=self.prefix[i]
        if i<len(self.starts) and self.starts[i]<x:
            answer += x-self.starts[i]
        return answer
    def overlap(self,start,end):return self.before(end)-self.before(start)

def fasta_ids(path):
    ids=set()
    for line in path.open():
        if line.startswith('>'):
            pid=line[1:].split()[0]
            if pid in ids:raise ValueError(f'Duplicate FASTA ID: {pid}')
            ids.add(pid)
    return ids

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def classify(signals, te_pct, unknown_pct):
    host=bool(signals.get('host_context'))
    if host and (signals.get('explicit_TE') or signals.get('RT_or_integrase')):
        return 'host_context_review'
    if signals.get('explicit_TE'):return 'priority_functional_review'
    if signals.get('RT_or_integrase'):return 'functional_candidate_review'
    if signals.get('generic_associated'):return 'generic_domain_review'
    if te_pct>=50 or unknown_pct>=50:return 'repeat_context_only'
    return 'no_screening_signal'

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/work/cyu/annotation/scotti-annotation'))
    p.add_argument('--assembly-root',type=Path,default=Path('/work/cyu/assembly/worm/scotti'))
    p.add_argument('--candidate',type=Path)
    p.add_argument('--repeat-out',type=Path)
    p.add_argument('--out',type=Path)
    p.add_argument('--expected-proteins',type=int,default=11464)
    p.add_argument('--expected-genes',type=int,default=9234)
    a=p.parse_args(argv)
    candidate=a.candidate or a.root/'results/filtered_candidates/Scotti_filtered7_v1'
    out=a.out or a.root/'results/annotation_qc/08_TE_protein_review_v1'
    if out.exists():raise FileExistsError(f'Output exists; inspect rather than overwrite: {out}')
    inputs={
      'protein':candidate/'Scotti_nuclear.filtered7.proteins.fa',
      'primary':candidate/'Scotti_nuclear.filtered7.longest.proteins.fa',
      'gff':candidate/'Scotti_nuclear.filtered7.candidate.gff3',
      'evidence':candidate/'functional_evidence_v1/Scotti_filtered7.all_isoforms.evidence.tsv',
      'interpro':candidate/'functional_evidence_v1/Scotti_filtered7.interproscan.tsv',
      'repeat_out':a.repeat_out or a.assembly_root/'repeats/02_repeatmasker/Scotti_nuclear.mt_filtered.ragtag.provisional.fa.out',
    }
    for name,path in inputs.items():
        if not path.is_file() or not path.stat().st_size:raise FileNotFoundError(f'{name}: {path}')
    keep=fasta_ids(inputs['protein']);primary=fasta_ids(inputs['primary'])
    assert len(keep)==a.expected_proteins, len(keep)
    assert len(primary)==a.expected_genes and primary<=keep
    with inputs['evidence'].open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'protein_id','gene_id','length_aa','swissprot_title','swissprot_evalue',
                  'swissprot_query_coverage_pct','eggnog_description','interpro_descriptions','eggnog_PFAMs'}
        assert required<=set(reader.fieldnames), required-set(reader.fieldnames)
        evidence={}
        for row in reader:
            pid=row['protein_id']
            if pid in evidence:raise ValueError(f'Duplicate evidence row: {pid}')
            evidence[pid]=row
    assert set(evidence)==keep
    genes={e['gene_id'] for e in evidence.values()}
    assert len(genes)==a.expected_genes
    assert {evidence[p]['gene_id'] for p in primary}==genes
    tx_gene={};cds=defaultdict(list);locations={}
    for line in inputs['gff'].open():
        if line.startswith('##FASTA'):break
        if not line.strip() or line.startswith('#'):continue
        f=line.rstrip('\r\n').split('\t')
        if len(f)!=9:raise ValueError('Malformed GFF3')
        # Split escaped values only after separating GFF delimiters.
        at=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        if f[2]=='mRNA':
            pid=unquote(at['ID']);gene=unquote(at['Parent'])
            if pid in tx_gene:raise ValueError('Duplicate transcript')
            tx_gene[pid]=gene;locations[pid]=(f[0],int(f[3]),int(f[4]),f[6])
        elif f[2]=='CDS':
            for parent in at.get('Parent','').split(','):
                pid=unquote(parent)
                cds[pid].append((f[0],int(f[3])-1,int(f[4])))
    assert set(tx_gene)==keep
    assert all(tx_gene[p]==evidence[p]['gene_id'] for p in keep)
    assert all(p in cds for p in keep)
    print('PASS: filtered7 protein/GFF/evidence IDs agree.',flush=True)

    signals={p:defaultdict(set) for p in keep}
    def scan(pid,source,value):
        if not value or value in {'-','none'}:return
        for name,regex in COMPILED.items():
            hits=sorted({m.group(0) for m in regex.finditer(value)})
            if hits:signals[pid][name].add(source+': '+', '.join(hits))
    for pid,e in evidence.items():
        for field in ['swissprot_title','eggnog_description','eggnog_PFAMs','interpro_descriptions']:
            scan(pid,field,e.get(field,''))
    ipr_match_rows=[]
    print('Scanning InterPro descriptions...',flush=True)
    for line in inputs['interpro'].open():
        if not line.strip() or line.startswith('#'):continue
        f=line.rstrip('\r\n').split('\t')
        if len(f)<11:raise ValueError('Malformed InterProScan row')
        pid=f[0]
        if pid not in keep:raise ValueError(f'Unknown InterPro query: {pid}')
        if int(f[2])!=int(evidence[pid]['length_aa']):
            raise ValueError(f'InterPro protein length mismatch: {pid}')
        descriptions=' | '.join([f[4],f[5]]+([f[11],f[12]] if len(f)>12 else []))
        matched=[k for k,r in COMPILED.items() if r.search(descriptions)]
        if matched:
            scan(pid,f'{f[3]}:{f[4]}',descriptions)
            ipr_match_rows.append([pid,tx_gene[pid],';'.join(matched),f[3],f[4],f[5],f[6],f[7],f[8],f[11] if len(f)>11 else '',f[12] if len(f)>12 else ''])

    print('Indexing RepeatMasker intervals...',flush=True)
    repeat_data=defaultdict(lambda:defaultdict(list));rm_rows=0
    for line in inputs['repeat_out'].open():
        f=line.split()
        if len(f)<11 or not f[0].isdigit():continue
        seq=f[4];s,e=int(f[5])-1,int(f[6]);cls=f[10]
        if s<0 or e<=s:raise ValueError('Invalid RepeatMasker interval')
        repeat_data[seq]['all'].append((s,e));rm_rows+=1
        root_class=cls.split('/')[0]
        if root_class in TE_CLASSES:repeat_data[seq]['classified_TE'].append((s,e))
        elif root_class.lower() in {'unknown','unclassified'}:repeat_data[seq]['unknown'].append((s,e))
    if not rm_rows:raise ValueError('No RepeatMasker records parsed')
    relevant={loc[0] for loc in locations.values()}
    missing=relevant-set(repeat_data)
    # A contig without a reported repeat is valid; record it explicitly.
    indexes={s:{k:Coverage(v) for k,v in groups.items()} for s,groups in repeat_data.items() if s in relevant}
    del repeat_data
    rows=[];by_gene=defaultdict(list)
    for pid in sorted(keep):
        e=evidence[pid];seqid,start,end,strand=locations[pid]
        segments=[]
        for seq,s,t in cds[pid]:
            if seq!=seqid:raise ValueError('Cross-contig CDS')
            segments.append((s,t))
        segments=merge(segments);bp=sum(t-s for s,t in segments)
        covered={}
        for category in ['all','classified_TE','unknown']:
            idx=indexes.get(seqid,{}).get(category)
            covered[category]=sum(idx.overlap(s,t) for s,t in segments) if idx else 0
        assert 0<=covered['classified_TE']<=covered['all']<=bp
        assert 0<=covered['unknown']<=covered['all']
        percent={k:100*v/bp for k,v in covered.items()}
        sig=signals[pid]
        tier=classify(sig,percent['classified_TE'],percent['unknown'])
        reasons=[]
        if percent['classified_TE']>=50:reasons.append('classified_TE_CDS_overlap_ge_50pct_context_only')
        if percent['unknown']>=50:reasons.append('unknown_repeat_CDS_overlap_ge_50pct_context_only')
        for k in PATTERNS:
            if sig.get(k):reasons.append(k)
        row=dict(protein_id=pid,gene_id=tx_gene[pid],primary_isoform=pid in primary,
                 length_aa=e['length_aa'],seqid=seqid,start=start,end=end,strand=strand,
                 CDS_union_bp=bp,CDS_RM_all_pct=round(percent['all'],4),
                 CDS_RM_classified_TE_pct=round(percent['classified_TE'],4),
                 CDS_RM_unknown_pct=round(percent['unknown'],4),review_tier=tier,
                 screening_reasons=';'.join(reasons),
                 explicit_TE_signals=' | '.join(sorted(sig.get('explicit_TE',set()))),
                 RT_integrase_signals=' | '.join(sorted(sig.get('RT_or_integrase',set()))),
                 generic_signals=' | '.join(sorted(sig.get('generic_associated',set()))),
                 host_context_signals=' | '.join(sorted(sig.get('host_context',set()))),
                 swissprot_title=e.get('swissprot_title',''),swissprot_evalue=e.get('swissprot_evalue',''),
                 swissprot_query_coverage_pct=e.get('swissprot_query_coverage_pct',''),
                 eggnog_description=e.get('eggnog_description',''),eggnog_PFAMs=e.get('eggnog_PFAMs',''),
                 interpro_descriptions=e.get('interpro_descriptions',''),
                 manual_decision='pending',reviewer_reason='')
        rows.append(row);by_gene[tx_gene[pid]].append(row)
    out.mkdir(parents=True)
    def write_table(name,items,fields=None):
        if fields is None:fields=list(rows[0])
        with (out/name).open('x') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(items)
    write_table('Scotti_TE_review.all_isoforms.tsv',rows)
    write_table('Scotti_TE_review.primary.tsv',[r for r in rows if r['primary_isoform']])
    write_table('Scotti_TE_review.flagged.tsv',[r for r in rows if r['review_tier']!='no_screening_signal'])
    functional_tiers={'priority_functional_review','functional_candidate_review','host_context_review'}
    shortlist=[r['protein_id'] for r in rows if r['primary_isoform'] and r['review_tier'] in functional_tiers]
    (out/'Scotti.primary.TE_functional_candidate_ids.txt').write_text(''.join(p+'\n' for p in shortlist))
    gene_rows=[]
    for gene,rs in sorted(by_gene.items()):
        pr=[r for r in rs if r['primary_isoform']];assert len(pr)==1
        gene_rows.append(dict(gene_id=gene,primary_protein=pr[0]['protein_id'],primary_tier=pr[0]['review_tier'],
            isoforms=len(rs),flagged_isoforms=';'.join(r['protein_id'] for r in rs if r['review_tier']!='no_screening_signal'),
            functional_candidate_isoforms=';'.join(r['protein_id'] for r in rs if r['review_tier'] in functional_tiers),
            any_host_context=any(bool(signals[r['protein_id']].get('host_context')) for r in rs),
            manual_decision='pending',reviewer_reason=''))
    write_table('Scotti_TE_review.genes.tsv',gene_rows,list(gene_rows[0]))
    with (out/'matched_interpro_evidence.tsv').open('x') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['protein_id','gene_id','signal_categories','database','signature','signature_description','aa_start','aa_end','score','interpro_id','interpro_description']);w.writerows(ipr_match_rows)
    with (out/'screening_rules.tsv').open('x') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['category','regex']);w.writerows(PATTERNS.items())
    summary=dict(status='candidate_review_only_no_sequences_removed',proteins=len(keep),genes=len(genes),primary_proteins=len(primary),
        all_isoform_tiers=dict(Counter(r['review_tier'] for r in rows)),
        primary_tiers=dict(Counter(r['review_tier'] for r in rows if r['primary_isoform'])),
        primary_functional_candidates=len(shortlist),repeatmasker_records=rm_rows,
        gene_contigs_without_reported_RM_hits=sorted(missing))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'README.txt').write_text(
      'Scotti filtered7 TE-protein candidate review; not a final TE classification.\n'
      'No proteins, genes, CDS or GFF features were removed or edited. No nonTE proteome was generated.\n'
      'Priority/functional candidate tiers use description/signature screening, not independent confirmation. Swiss-Prot and eggNOG hits can propagate similar homology evidence. Inspect matched_interpro_evidence.tsv for individual domains, positions and scores.\n'
      'RT/integrase and RNase H can occur in host proteins. Host context is flagged and requires manual review, not automatic retention or exclusion. Generic nuclease/protease signals are weak and never confirm TE origin.\n'
      'RepeatMasker overlaps are union lengths on CDS, not gene spans or a direct measurement of lowercase FASTA. Classified TE overlap excludes Unknown, simple repeats, low complexity, satellites and small RNA. Unknown repeats are reported separately.\n'
      'The 50% overlap threshold is only a review-priority heuristic; overlap does not establish TE derivation. Repeat-only and no-signal models must not be automatically labelled nonTE. Small RNA was not a TE class.\n'
      'Candidate ID list contains primary proteins with functional screening signals, including host-context review cases. It is NOT an exclusion list. Gene table flags signal in alternative isoforms too. Manual decisions start as pending.\n'
      'Source files are the existing filtered7 evidence package and the matching RagTag RepeatMasker output.\n')
    script=Path(__file__).resolve()
    shutil.copy2(script,out/'build_scotti_te_review.py')
    print('Recording input checksums...',flush=True)
    with (out/'input.SHA256SUMS').open('x') as f:
        for path in inputs.values():f.write(f'{digest(path)}  {path}\n')
    with (out/'SHA256SUMS').open('x') as f:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.name!='SHA256SUMS':f.write(f'{digest(path)}  {path.name}\n')
    print(json.dumps(summary,indent=2),flush=True);print('Output:',out,flush=True)

if __name__=='__main__':main()
