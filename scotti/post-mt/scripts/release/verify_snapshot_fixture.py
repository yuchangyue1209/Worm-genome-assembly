#!/usr/bin/env python3
"""Meaningful synthetic checks; real GenomeTools validation runs on the server."""
import tempfile
from pathlib import Path
import freeze_scotti_working_snapshot as m


def rejects(fn, text):
    try:
        fn()
    except ValueError as e:
        assert text in str(e), str(e)
    else:
        raise AssertionError('Expected rejection: '+text)


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    f = root/'multipart.gff3'
    f.write_text('##gff-version 3\n'
        'chr\tBRAKER\tgene\t1\t100\t.\t+\t.\tID=g1\n'
        'chr\tBRAKER\tmRNA\t1\t100\t.\t+\t.\tID=g1.t1;Parent=g1\n'
        'chr\tBRAKER\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1.t1\n'
        'chr\tBRAKER\tCDS\t51\t59\t.\t+\t0\tID=c1;Parent=g1.t1\n')
    assert len(m.gff(f,{'chr':100}))==4
    f.write_text(f.read_text().replace('Parent=g1.t1','Parent=absent'))
    rejects(lambda:m.gff(f,{'chr':100}), 'Unresolved')
    f.write_text('##gff-version 3\nchr\tx\tgene\t1\t101\t.\t+\t.\tID=g\n')
    rejects(lambda:m.gff(f,{'chr':100}), 'coordinate mismatch')
    f.write_text('>same\nMAA\n>same\nMAA\n')
    rejects(lambda:m.fasta(f), 'Duplicate FASTA ID')
    f.write_text('MAA\n')
    rejects(lambda:m.fasta(f), 'Sequence before FASTA')
    full = root/'full_table.tsv'
    rows = []
    counter = 0
    for status,n in [('Complete',487),('Duplicated',7),('Fragmented',66),('Missing',112)]:
        for i in range(n):
            counter += 1
            rows.append(f'{counter}at33208\t{status}\tp1\n')
            if status=='Duplicated':rows.append(f'{counter}at33208\t{status}\tp2\n')
    full.write_text('# BUSCO fixture\n'+''.join(rows))
    states = m.busco_states(full)
    assert len(states)==672 and m.Counter(states.values())['Duplicated']==7

# Verify the production count and one-per-gene logic with realistic set sizes.
features=[]; all_proteins={}; primary={}; coding={}
for i in range(9234):
    gene=f'gene{i:05d}'; pid=gene+'.t1'
    features.extend([(['chr','BRAKER','gene','1','9','.','+','.',''],{'ID':gene}),
                     (['chr','BRAKER','mRNA','1','9','.','+','.',''],{'ID':pid,'Parent':gene})])
    all_proteins[pid]='MAAA'; primary[pid]='MAAA'; coding[pid]='ATGGCTGCTGCT'
    if i<2230:
        alt=gene+'.t2'
        features.append((['chr','BRAKER','mRNA','1','9','.','+','.',''],{'ID':alt,'Parent':gene}))
        all_proteins[alt]='MA'; coding[alt]='ATGGCT'
genes, tx=m.check_proteins(features,all_proteins,coding,primary)
assert len(genes)==9234 and len(tx)==11464
broken=dict(primary);broken.pop(next(iter(broken)))
rejects(lambda:m.check_proteins(features,all_proteins,coding,broken),'primary protein set')
broken=dict(primary);p=next(iter(broken));broken[p]='MMMM'
rejects(lambda:m.check_proteins(features,all_proteins,coding,broken),'sequence mismatch')
broken=dict(primary);p='gene00000.t1';broken.pop(p);broken['gene00000.t2']='MA'
rejects(lambda:m.check_proteins(features,all_proteins,coding,broken),'not a longest isoform')
print('PASS: GFF3 coordinates/parents and multipart CDS IDs; FASTA duplicate detection;')
print('      production gene/isoform counts, longest representatives and BUSCO duplicate grouping.')
