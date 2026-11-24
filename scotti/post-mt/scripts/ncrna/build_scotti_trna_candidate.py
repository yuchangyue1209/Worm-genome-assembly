#!/usr/bin/env python3
"""Export high-confidence tRNA candidates and audit CDS overlaps.
No source files or protein models are edited. Python standard library only.
"""
from pathlib import Path
from collections import Counter, defaultdict
from bisect import bisect_right
from urllib.parse import quote
import argparse, csv, hashlib, json, re, shutil

class Coverage:
    def __init__(self, intervals):
        merged=[]
        for s,e in sorted(intervals):
            if merged and s<=merged[-1][1]:merged[-1]=(merged[-1][0],max(e,merged[-1][1]))
            else:merged.append((s,e))
        self.starts=[s for s,e in merged];self.ends=[e for s,e in merged];self.prefix=[0]
        for s,e in merged:self.prefix.append(self.prefix[-1]+e-s)
    def before(self,x):
        i=bisect_right(self.ends,x);n=self.prefix[i]
        if i<len(self.starts) and self.starts[i]<x:n+=x-self.starts[i]
        return n
    def overlap(self,s,e):return self.before(e)-self.before(s)

def lengths(path):
    result={};name=None
    for line in path.open():
        if line.startswith('>'):
            name=line[1:].split()[0]
            if name in result:raise ValueError(f'Duplicate FASTA ID: {name}')
            result[name]=0
        elif line.strip():
            if name is None:raise ValueError('Sequence without header')
            result[name]+=len(''.join(line.split()))
    return result

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def attrs(values):
    return ';'.join(k+'='+quote(str(v),safe='._:-') for k,v in values.items())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/work/cyu/annotation/scotti-annotation'))
    p.add_argument('--genome',type=Path,default=Path('/work/cyu/assembly/worm/scotti/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa'))
    p.add_argument('--out',type=Path)
    p.add_argument('--expected',type=int,default=344)
    a=p.parse_args()
    base=a.root/'results/ncrna_v1'
    table=base/'03b_trna_highconfidence_breakdown/Scotti_tRNA.out'
    log=base/'03b_trna_highconfidence_breakdown/Scotti_tRNA.log'
    pc=a.root/'results/filtered_candidates/Scotti_filtered7_v1/Scotti_nuclear.filtered7.candidate.gff3'
    rfam=base/'02_rfam/Scotti_nuclear.rfam.tblout'
    out=a.out or base/'04_trna_candidate_v1'
    if out.exists():raise FileExistsError(f'Output exists; inspect before rerunning: {out}')
    for path in [table,log,pc,rfam,a.genome]:
        if not path.is_file() or not path.stat().st_size:raise FileNotFoundError(path)
    with table.open() as f:header=[next(f) for _ in range(3)]
    header_text=''.join(header)
    assert 'HMM' in header_text and "2'Str" in header_text, 'Missing score breakdown; do not filter without these columns.'
    name_header=next((l for l in header if l.lstrip().startswith('Name')),None)
    assert name_header is not None
    cols=[v.strip() for v in name_header.replace('tRNA #','tRNA#').rstrip('\r\n').split('\t')]
    assert 'Note' in cols and 'CM' in cols
    note_i=cols.index('Note');cm_i=cols.index('CM')
    iso_score_i=cm_i+1
    score_headers=[v.strip() for v in header[0].replace('Intron Bounds','Intron\tBound').rstrip('\r\n').split('\t')]
    assert 'Inf' in score_headers
    score_i=score_headers.index('Inf')
    hmm_i=score_headers.index('HMM');ss_i=score_headers.index("2'Str")
    # Ensure the breakdown is situated where both table readers expect it.
    assert score_i<cm_i and hmm_i<cm_i and ss_i<cm_i
    seq_lengths=lengths(a.genome)
    rows=[];mismatches=[];total=0;prediction_ids=set()
    with table.open() as f:
        for line in f:
            flds=[v.strip() for v in line.rstrip('\r\n').split('\t')]
            if len(flds)<2 or not flds[1].isdigit():continue
            if len(flds)<=max(note_i,iso_score_i,score_i,ss_i,hmm_i):raise ValueError(f'Short row: {line!r}')
            total+=1;pid=(flds[0],int(flds[1]))
            if pid in prediction_ids:raise ValueError(f'Duplicate prediction ID: {pid}')
            prediction_ids.add(pid)
            tags={t.strip() for t in flds[note_i].split(',') if t.strip()}
            begin,end=int(flds[2]),int(flds[3]);s,e=sorted([begin,end])
            if flds[0] not in seq_lengths or not 1<=s<=e<=seq_lengths[flds[0]]:
                raise ValueError(f'Invalid genome coordinates: {pid}')
            r=dict(seqid=flds[0],prediction_number=flds[1],start=s,end=e,
                strand='+' if begin<=end else '-',isotype=flds[4],anticodon=flds[5],
                infernal_score=float(flds[score_i]),HMM_score=float(flds[hmm_i]),
                structure_score=float(flds[ss_i]),isotype_model=flds[cm_i],
                isotype_score=float(flds[iso_score_i]),original_notes=flds[note_i],
                intron_columns_original=' | '.join(flds[6:score_i]))
            if 'high confidence set' in tags:
                assert 'pseudo' not in tags and 'isotype mismatch' not in tags
                rows.append(r)
            elif 'isotype mismatch' in tags:mismatches.append(r)
    expected_match=re.search(r'^High confidence set:\s*(\d+)',log.read_text(),re.M)
    assert expected_match, 'High-confidence count missing from official filter log'
    assert len(rows)==int(expected_match.group(1))==a.expected, (len(rows),expected_match.group(1),a.expected)
    rows.sort(key=lambda r:(r['seqid'],r['start'],r['end'],r['strand']))
    unique={(r['seqid'],r['start'],r['end'],r['strand']) for r in rows}
    assert len(unique)==len(rows), 'Exact duplicate loci require manual review'

    cds=defaultdict(list);cds_count=0
    for line in pc.open():
        if line.startswith('##FASTA'):break
        if not line.strip() or line.startswith('#'):continue
        f=line.rstrip().split('\t')
        if len(f)!=9:raise ValueError('Malformed protein GFF3')
        if f[2]!='CDS':continue
        s,e=int(f[3])-1,int(f[4])
        assert f[0] in seq_lengths and 0<=s<e<=seq_lengths[f[0]]
        cds[(f[0],'all')].append((s,e));cds[(f[0],f[6])].append((s,e));cds_count+=1
    assert cds_count>0
    indexes={key:Coverage(intervals) for key,intervals in cds.items()}
    for i,r in enumerate(rows,1):
        r['feature_id']=f'Scotti_tRNA_{i:04d}'
        for strand,label in [('all','CDS_overlap_bp'),(r['strand'],'same_strand_CDS_overlap_bp')]:
            index=indexes.get((r['seqid'],strand))
            r[label]=index.overlap(r['start']-1,r['end']) if index else 0
        r['review_note']='predicted_SeC_function_not_validated' if r['isotype'] in {'SeC','SelCys','Sec'} else 'model_predicted_high_confidence'
    # Audit Rfam families without selecting/merging ncRNA models yet.
    families=defaultdict(Counter);rfam_n=0
    for line in rfam.open():
        if not line.strip() or line.startswith('#'):continue
        f=line.split(maxsplit=28)
        assert len(f)>=28 and f[0].isdigit() and f[2].startswith('RF'), 'Unexpected Rfam fmt2 record'
        assert f[3] in seq_lengths and int(f[27])==seq_lengths[f[3]], 'Rfam query/genome mismatch'
        assert f[11] in {'+','-'} and 1<=min(int(f[9]),int(f[10]))<=max(int(f[9]),int(f[10]))<=seq_lengths[f[3]]
        rfam_n+=1;k=(f[1],f[2]);families[k]['hits']+=1
        if f[18]=='!':families[k]['included_hits']+=1
        if f[12]=='no':families[k]['untruncated_hits']+=1
        families[k]['olp_'+f[19]]+=1
    assert rfam_n>0
    out.mkdir(parents=True)
    with (out/'Scotti_nuclear.tRNA.high_confidence.candidate.gff3').open('x') as g, (out/'Scotti_nuclear.tRNA.high_confidence.candidate.bed').open('x') as b:
        g.write('##gff-version 3\n')
        for seq in sorted({r['seqid'] for r in rows}):g.write(f'##sequence-region {seq} 1 {seq_lengths[seq]}\n')
        for r in rows:
            attributes=attrs(dict(ID=r['feature_id'],Name='tRNA-'+r['isotype']+'-'+r['anticodon'],
                product='tRNA-'+r['isotype'],isotype=r['isotype'],anticodon=r['anticodon'],
                prediction_number=r['prediction_number'],isotype_score=r['isotype_score'],
                structure_score=r['structure_score'],Note='high confidence set; '+r['review_note']))
            g.write('\t'.join(map(str,[r['seqid'],'tRNAscan-SE','tRNA',r['start'],r['end'],r['infernal_score'],r['strand'],'.',attributes]))+'\n')
            b.write('\t'.join(map(str,[r['seqid'],r['start']-1,r['end'],r['feature_id'],'.',r['strand']]))+'\n')
    fields=list(rows[0])
    with (out/'tRNA_high_confidence.audit.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    with (out/'tRNA_isotype_mismatch.review.tsv').open('x') as f:
        mismatch_fields=[k for k in fields if k not in {'feature_id','CDS_overlap_bp','same_strand_CDS_overlap_bp','review_note'}]
        w=csv.DictWriter(f,fieldnames=mismatch_fields,delimiter='\t');w.writeheader();w.writerows(mismatches)
    with (out/'Rfam_family_summary.tsv').open('x') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['family','accession','hits','included_hits','untruncated_hits','overlap_annotation_counts'])
        for (family,accession),c in sorted(families.items(),key=lambda item:(-item[1]['hits'],item[0])):
            w.writerow([family,accession,c['hits'],c['included_hits'],c['untruncated_hits'],json.dumps({k[4:]:v for k,v in c.items() if k.startswith('olp_')},sort_keys=True)])
    summary=dict(status='tRNA_candidate_only_not_merged_into_protein_annotation',
        source_tagged_predictions=total,high_confidence_tRNA_features=len(rows),
        isotype_mismatch_review=len(mismatches),isotypes=dict(sorted(Counter(r['isotype'] for r in rows).items())),
        tRNA_features_overlapping_any_CDS=sum(r['CDS_overlap_bp']>0 for r in rows),
        tRNA_features_overlapping_same_strand_CDS=sum(r['same_strand_CDS_overlap_bp']>0 for r in rows),
        Rfam_raw_hits=rfam_n,Rfam_families=len(families),
        source_protein_models_removed=0)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'README.txt').write_text('Candidate nuclear tRNA annotation derived ONLY from high confidence set tags in the corrected score-breakdown tRNAscan-SE/EukHighConfidenceFilter output. Official log count and FASTA coordinates validated. One standalone tRNA feature per precursor locus, consistent with the earlier flat ncRNA workflow; exon/intron features and mature tRNA sequences are not constructed here. Intron information is retained verbatim in audit.tsv. No genome, protein models or source predictions were changed. CDS overlap is reported as union-covered bases on the tRNA genomic span, across all strands and same strand; overlap alone is not a removal criterion. Isotype-mismatch models are exported separately and not merged into the high-confidence candidate. SeC is a computational prediction with unverified function. Rfam family counts are raw-hit diagnostics, not final gene counts. Rfam has not been deduplicated or merged into this GFF3. External GFF validation is a separate step.\n')
    shutil.copy2(Path(__file__),out/'build_scotti_trna_candidate.py')
    with (out/'input.SHA256SUMS').open('x') as f:
        for path in [table,log,pc,rfam,a.genome]:f.write(f'{sha(path)}  {path}\n')
    with (out/'SHA256SUMS').open('x') as f:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.name!='SHA256SUMS':f.write(f'{sha(path)}  {path.name}\n')
    print(json.dumps(summary,indent=2));print('Output:',out)

if __name__=='__main__':main()
