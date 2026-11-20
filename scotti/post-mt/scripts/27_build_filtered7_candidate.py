#!/usr/bin/env python3
"""Build a reversible candidate excluding seven reviewed extremely short models.
Requires Biopython. Does not edit source predictions or merge neighboring genes.
"""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import json
from Bio import SeqIO

REASONS = {
    'g1218': '10 aa; no functional hits; fully softmasked CDS; RNA signal does not validate ORF',
    'g196': '5 aa; no functional hits; 88 bp from same-strand g195; retain locus for repair review',
    'g3336': '8 aa; no functional hits; 22-bp CDS; 61 bp from same-strand g3337; retain locus for repair review',
    'g5529': '10 aa; no functional hits; 210 bp from same-strand g5528; strong regional RNA; retain locus for repair review',
    'g5573': '4 aa; no functional hits; fully softmasked CDS; no unique regional RNA records',
    'g7051': '4 aa; no functional hits; gene span contained in same-strand g7050; no automatic merging',
    'g7437': '2 aa; no functional hits; 9-bp CDS; sparse regional RNA records',
}

def attrs(text):
    return dict(item.split('=', 1) for item in text.split(';') if '=' in item)

def load_fasta(path):
    result = {}
    for record in SeqIO.parse(path, 'fasta'):
        if record.id in result:
            raise ValueError(f'Duplicate FASTA ID: {record.id}')
        result[record.id] = record
    return result

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, default=Path('/work/cyu/annotation/scotti-annotation/results/braker3_etp_v1_retry2'))
    p.add_argument('--out', type=Path, default=Path('/work/cyu/annotation/scotti-annotation/results/filtered_candidates/Scotti_filtered7_v1'))
    p.add_argument('--expected-genes', type=int, default=9241)
    p.add_argument('--expected-transcripts', type=int, default=11471)
    a = p.parse_args()
    if a.out.exists():
        raise FileExistsError(f'Output already exists: {a.out}')
    raw = (a.run/'braker.gff3').read_text().splitlines(keepends=True)
    features = []
    genes, tx_gene = {}, {}
    ids = set()
    for line in raw:
        if line.startswith('##FASTA'):
            raise ValueError('Embedded FASTA requires separate handling')
        if not line.strip() or line.startswith('#'):
            features.append((line, None, None)); continue
        f = line.rstrip('\r\n').split('\t')
        if len(f) != 9:
            raise ValueError('Malformed GFF3 record')
        at = attrs(f[8]); features.append((line, f, at))
        if at.get('ID'):
            ids.add(at['ID'])
        if f[2] == 'gene':
            if at['ID'] in genes: raise ValueError('Duplicate gene ID')
            genes[at['ID']] = f
        elif f[2] == 'mRNA':
            if at['ID'] in tx_gene or ',' in at['Parent']: raise ValueError('Unexpected mRNA parents')
            tx_gene[at['ID']] = at['Parent']
    assert len(genes) == a.expected_genes, len(genes)
    assert len(tx_gene) == a.expected_transcripts, len(tx_gene)
    assert set(tx_gene.values()) == set(genes), 'Gene/parent mismatch'
    assert set(REASONS) <= genes.keys(), 'Reviewed genes missing'
    proteins = load_fasta(a.run/'braker.aa')
    coding = load_fasta(a.run/'braker.codingseq')
    assert set(proteins) == set(coding) == set(tx_gene), 'GFF/protein/CDS ID mismatch'
    removed_tx = {tx for tx, gene in tx_gene.items() if gene in REASONS}
    assert len(removed_tx) == 7, 'Unexpected isoform count in reviewed genes'
    assert all(len(str(proteins[tx].seq).rstrip('*')) <= 10 for tx in removed_tx)
    removed_ids = set(REASONS) | removed_tx
    changed = True
    while changed:
        changed = False
        for _, f, at in features:
            if f is None: continue
            parents = set(filter(None, at.get('Parent', '').split(',')))
            if parents & removed_ids:
                if parents - removed_ids:
                    raise ValueError('Mixed retained/excluded parent feature needs manual handling')
                if at.get('ID') and at['ID'] not in removed_ids:
                    removed_ids.add(at['ID']); changed = True
    kept = []
    counts = Counter()
    for line, f, at in features:
        if f is None:
            # Feature-group terminators are dispensable after filtering.
            if line.strip() != '###': kept.append(line)
            continue
        parents = set(filter(None, at.get('Parent', '').split(',')))
        if at.get('ID') in removed_ids or parents & removed_ids: continue
        if parents - ids: raise ValueError('Unknown GFF parent')
        kept.append(line); counts[f[2]] += 1
    assert counts['gene'] == a.expected_genes-7
    assert counts['mRNA'] == a.expected_transcripts-7
    remaining = [r for tx, r in proteins.items() if tx not in removed_tx]
    remaining_cds = [r for tx, r in coding.items() if tx not in removed_tx]
    by_gene = defaultdict(list)
    for r in remaining: by_gene[tx_gene[r.id]].append(r)
    longest = [min(rs, key=lambda r: (-len(str(r.seq).rstrip('*')), r.id)) for _, rs in sorted(by_gene.items())]
    a.out.mkdir(parents=True)
    # GenomeTools requires the version directive on the very first line.
    body = []
    for line in kept:
        line = line.lstrip('\ufeff')
        if line.startswith('##gff-version'):
            if line.strip() != '##gff-version 3':
                raise ValueError('Unexpected GFF version directive')
            continue
        body.append(line)
    (a.out/'Scotti_nuclear.filtered7.candidate.gff3').write_text(
        '##gff-version 3\n' + ''.join(body)
    )
    SeqIO.write(remaining, a.out/'Scotti_nuclear.filtered7.proteins.fa', 'fasta')
    SeqIO.write(remaining_cds, a.out/'Scotti_nuclear.filtered7.cds.fa', 'fasta')
    SeqIO.write(longest, a.out/'Scotti_nuclear.filtered7.longest.proteins.fa', 'fasta')
    SeqIO.write([proteins[tx] for tx in sorted(removed_tx)], a.out/'excluded.proteins.fa', 'fasta')
    with (a.out/'excluded_models.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['gene_id','transcript_id','protein_aa','seqid','start','end','strand','reason'])
        for tx in sorted(removed_tx):
            gene=tx_gene[tx];g=genes[gene]
            w.writerow([gene,tx,len(str(proteins[tx].seq).rstrip('*')),g[0],g[3],g[4],g[6],REASONS[gene]])
    with (a.out/'longest_isoform_selection.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['gene_id','transcript_id'])
        for r in longest:w.writerow([tx_gene[r.id],r.id])
    stats=dict(genes=counts['gene'],transcripts=counts['mRNA'],proteins=len(remaining),cds_sequences=len(remaining_cds),longest_proteins=len(longest),excluded_genes=sorted(REASONS),status='candidate_not_final_release')
    (a.out/'summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    (a.out/'README.txt').write_text('Candidate only: excludes seven reviewed 2–10 aa protein-coding models.\nOriginal BRAKER files and genome sequence are unchanged. No neighboring genes were merged.\nExclusion does not establish biological gene absence. Retain loci for possible model repair, particularly g196, g3336 and g5529.\nThe remaining 37 short-gene candidates are retained; the two short alternative isoforms in longer genes are retained.\nNo genome coordinates were altered. Remaining proteins/CDS are copied from BRAKER without retranslation.\nExternal GFF3 validation and candidate BUSCO checks remain to be run.\n')
    with (a.out/'input.SHA256SUMS').open('w') as f:
        for path in [a.run/'braker.gff3', a.run/'braker.aa', a.run/'braker.codingseq']:
            f.write(f'{digest(path)}  {path}\n')
    with (a.out/'SHA256SUMS').open('w') as f:
        for path in sorted(a.out.iterdir()):
            if path.is_file() and path.name!='SHA256SUMS':f.write(f'{digest(path)}  {path.name}\n')
    print(json.dumps(stats,indent=2));print('Output:',a.out)

if __name__=='__main__': main()
