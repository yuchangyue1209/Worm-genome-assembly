#!/usr/bin/env python3
"""S. cotti mt workflow helpers; coordinate-specific audits are for the 15437-bp assembly.
Biopython is required only for extract-gff, extract-genbank and translation.
"""
import argparse
import gzip
import random
import re
from collections import Counter, defaultdict
from pathlib import Path


def fasta(path):
    records = {}
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f'Duplicate FASTA ID: {name}')
                records[name] = []
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                records[name].append(line.strip())
    return {k: ''.join(v) for k, v in records.items()}


def write_fasta(handle, name, seq):
    handle.write(f'>{name}\n')
    for i in range(0, len(seq), 80):
        handle.write(seq[i:i+80] + '\n')


def fastq(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    limit = 10_000_000
    with opener(path, 'rb') as handle:
        count = 0
        while True:
            lines = []
            for i in range(4):
                line = handle.readline(limit + 1)
                if i == 0 and not line:
                    return
                if not line or len(line) > limit:
                    raise ValueError(f'Record {count+1}: truncated record or line exceeds 10 MB')
                lines.append(line.rstrip(b'\r\n'))
            h, s, p, q = lines
            if not h.startswith(b'@') or not p.startswith(b'+') or not s or len(s) != len(q):
                raise ValueError(f'Invalid FASTQ record {count+1}')
            count += 1
            yield lines


def union_length(intervals):
    end = -1
    total = 0
    for start, stop in sorted(intervals):
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def main(a):
    if a.command == 'fastq-check':
        count = bases = longest = 0
        for h, s, p, q in fastq(a.input):
            count += 1
            bases += len(s)
            longest = max(longest, len(s))
            if count % 1_000_000 == 0:
                print(f'Checked {count:,} reads', flush=True)
        if not count or (a.expected is not None and count != a.expected):
            raise ValueError(f'Unexpected read count: {count}; expected {a.expected}')
        print(f'PASS\nReads: {count:,}\nTotal bases: {bases:,}\nMean length: {bases/count:.1f}\nLongest: {longest:,}')
    elif a.command == 'sample':
        records = list(fastq(a.input))
        selected = random.Random(a.seed).sample(records, min(a.number, len(records)))
        if not selected:
            raise ValueError('Empty input')
        with open(a.output, 'xb') as out:
            for record in selected:
                out.write(b'\n'.join(record) + b'\n')
        bases = sum(len(r[1]) for r in selected)
        print(f'Seed: {a.seed}; reads: {len(selected)}; bases: {bases}; nominal depth at 16kb: {bases/16000:.1f}')
    elif a.command == 'summarize':
        groups = defaultdict(list)
        with open(a.input) as f:
            for line in f:
                p = line.split()
                groups[(p[1], p[0])].append(p)
        print('scaffold\tlength_bp\tmt_ref\tHSPs\tmt_ref_covered_bp\tmt_ref_covered_pct\tscaffold_covered_bp\tscaffold_covered_pct\tfirst_hit\tlast_hit')
        for (sid, qid), rows in groups.items():
            qlen, slen = map(int, rows[0][10:12])
            qi = [(min(int(p[4]), int(p[5]))-1, max(int(p[4]), int(p[5]))) for p in rows]
            si = [(min(int(p[6]), int(p[7]))-1, max(int(p[6]), int(p[7]))) for p in rows]
            qc, sc = union_length(qi), union_length(si)
            print(sid, slen, qid, len(rows), qc, round(qc/qlen*100, 2), sc, round(sc/slen*100, 2), min(x[0] for x in si)+1, max(x[1] for x in si), sep='\t')
    elif a.command == 'split-candidates':
        seqs = fasta(a.input)
        ids = {'scaffold206', 'scaffold312', 'scaffold349', 'scaffold410'}
        if not ids <= seqs.keys():
            raise ValueError(f'Missing expected IDs: {ids-seqs.keys()}')
        out = Path(a.output)
        with open(out/'scotti_mt_candidates.fa', 'x') as mt, open(out/'Scotti_nuclear.mt_filtered.provisional.fa', 'x') as nuc:
            for name, seq in seqs.items():
                write_fasta(mt if name in ids else nuc, name, seq)
        total = sum(map(len, seqs.values()))
        mtlen = sum(len(seqs[k]) for k in ids)
        print(f'Input: {len(seqs)} sequences, {total} bp; nuclear: {len(seqs)-4}, {total-mtlen} bp; candidates: 4, {mtlen} bp')
        for name in sorted(ids):
            seq = seqs[name].upper()
            runs = [f'{m.start()+1}-{m.end()}' for m in re.finditer('N{10,}', seq)]
            print(name, len(seq), seq.count('N'), ','.join(runs) or 'none', sep='\t')
    elif a.command == 'rotate':
        records = fasta(a.input)
        if len(records) != 1:
            raise ValueError('Expected one mt contig')
        seq = next(iter(records.values()))
        cut = len(seq)//2
        with open(a.output, 'x') as out:
            write_fasta(out, 'mt_rotated', seq[cut:]+seq[:cut])
        print(f'Length: {len(seq)}; original junction boundary (0-based): {len(seq)-cut}')
    elif a.command in ('junction', 'repeat'):
        spanning, indels = set(), set()
        sizes = Counter()
        with open(a.input) as f:
            for line in f:
                p = line.rstrip().split('\t')
                length, start, end = int(p[6]), int(p[7]), int(p[8])
                if length != 15437:
                    raise ValueError('This audit window is specific to the 15437-bp assembly')
                if int(p[10]) == 0 or int(p[9])/int(p[10]) < .95 or int(p[11]) < 20:
                    continue
                if a.command == 'junction':
                    j = length-length//2
                    if start <= j-1000 and end >= j+1000:
                        spanning.add(p[0])
                    continue
                if not (start <= 7300 and end >= 11300):
                    continue
                cigar = next((x[5:] for x in p[12:] if x.startswith('cg:Z:')), None)
                if cigar is None:
                    continue
                spanning.add(p[0])
                pos = start
                for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
                    n = int(n)
                    if op == 'I':
                        if 8300 <= pos <= 10300 and n >= 100:
                            indels.add(p[0]); sizes[('insertion', n)] += 1
                    elif op in 'MDN=X':
                        if op == 'D' and pos < 10300 and pos+n > 8300 and n >= 100:
                            indels.add(p[0]); sizes[('deletion', n)] += 1
                        pos += n
        print(f'Qualifying spanning reads: {len(spanning)}')
        if a.command == 'repeat':
            print(f'Reads with >=100bp indel: {len(indels)}; without: {len(spanning-indels)}')
            for (kind, size), n in sizes.most_common(10):
                print(kind, size, n)
        print('Filters: alignment-wide identity >=95%, MAPQ >=20. Conditional counts, not a heteroplasmy estimate.')
    elif a.command == 'extract-gff':
        from Bio.Seq import Seq
        seqs = fasta(a.fasta)
        wanted = set(a.genes)
        hits = set()
        with open(a.gff) as f:
            for line in f:
                if line.startswith('##FASTA'):
                    break
                if line.startswith('#'):
                    continue
                p = line.rstrip().split('\t')
                if len(p) < 9 or p[2] not in ('gene', 'ncRNA_gene', 'rRNA', 'tRNA'):
                    continue
                attrs = dict(x.split('=', 1) for x in p[8].split(';') if '=' in x)
                name = attrs.get('Name', '').lower().split('(')[0]
                if name in wanted:
                    hits.add((name, p[0], int(p[3]), int(p[4]), p[6]))
        if {x[0] for x in hits} != wanted:
            raise ValueError(f'Not all requested genes found: {wanted}')
        with open(a.output, 'w') as out:
            for name, chrom, start, end, strand in sorted(hits):
                seq = Seq(seqs[chrom][start-1:end])
                if strand == '-':
                    seq = seq.reverse_complement()
                write_fasta(out, f'solidus_{name}_{start}_{end}', str(seq))
                print(name, start, end, strand, len(seq))
    elif a.command == 'extract-genbank':
        from Bio import SeqIO
        rec = SeqIO.read(a.input, 'genbank')
        if rec.id != 'MW602521.1' or rec.annotations.get('organism') != 'Schistocephalus solidus':
            raise ValueError('Unexpected accession/species')
        print(rec.id, rec.description)
        found = set()
        with open(a.output, 'w') as out:
            for feature in rec.features:
                gene = ' '.join(feature.qualifiers.get('gene', [])).lower()
                product = ' '.join(feature.qualifiers.get('product', [])).lower()
                if feature.type == 'tRNA' and (gene in ('trnh', 'trn-h') or 'his' in product):
                    label = 'trnH'
                elif feature.type == 'CDS' and (gene in ('cox3', 'co3', 'coiii') or 'oxidase subunit iii' in product or 'oxidase subunit 3' in product):
                    label = 'cox3'
                else:
                    continue
                seq = feature.extract(rec.seq)
                print(label, feature.location, len(seq), feature.qualifiers)
                write_fasta(out, f'{rec.id}_{label}', str(seq))
                found.add(label)
        if found != {'trnH', 'cox3'}:
            raise ValueError('Expected trnH and cox3 features')
    elif a.command == 'translation':
        from Bio.Seq import Seq
        cot = next(iter(fasta(a.input).values()))
        sol = next(iter(fasta(a.solidus).values()))
        for label, seq, start, end in [
            ('solidus_original', sol, 10469, 11114),
            ('cotti_MITOS2', cot, 5366, 6016),
            ('cotti_own_solidus_BLAST', cot, 5372, 6016),
            ('cotti_published_homology_candidate', cot, 5428, 6070),
        ]:
            dna = str(Seq(seq[start-1:end]).reverse_complement()).upper()
            n = len(dna)//3*3
            protein = str(Seq(dna[:n]).translate(table=9))
            print(f'\n{label}: {start}-{end}, minus, {len(dna)} bp')
            print('First30:', dna[:30], 'Last30:', dna[-30:], 'Remainder:', dna[n:] or 'none')
            print('Stop positions:', [i+1 for i, aa in enumerate(protein) if aa == '*'])
            print(protein)
        print('Direct translation is diagnostic: initial GTG prints V. Do not add genomic As to complete terminal T.')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for cmd in ('fastq-check', 'sample', 'summarize', 'split-candidates', 'rotate', 'junction', 'repeat', 'extract-genbank', 'translation'):
        s = sub.add_parser(cmd)
        s.add_argument('input')
        if cmd in ('sample', 'split-candidates', 'rotate', 'extract-genbank'):
            s.add_argument('output')
        if cmd == 'fastq-check':
            s.add_argument('--expected', type=int)
        if cmd == 'sample':
            s.add_argument('--number', type=int, default=700)
            s.add_argument('--seed', type=int, default=42)
        if cmd == 'translation':
            s.add_argument('solidus')
    s = sub.add_parser('extract-gff')
    s.add_argument('fasta'); s.add_argument('gff'); s.add_argument('output')
    s.add_argument('--genes', nargs='+', required=True)
    return p


if __name__ == '__main__':
    main(parser().parse_args())
