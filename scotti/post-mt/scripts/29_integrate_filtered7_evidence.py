#!/usr/bin/env python3
"""Filter existing functional TSVs by candidate IDs and build an evidence matrix.
Standard library only. No search reruns, gene-model edits, or product-name assignment.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import hashlib
import json
import shutil


def fasta_lengths(path):
    result = {}; name = None; sequence = []
    def save():
        if name is None: return
        if name in result: raise ValueError(f'Duplicate FASTA ID: {name}')
        seq = ''.join(sequence)
        if seq.endswith('*'): seq = seq[:-1]
        result[name] = len(seq)
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                save(); name = line[1:].split()[0]; sequence = []
            elif line.strip():
                if name is None: raise ValueError('Sequence before FASTA header')
                sequence.append(line.strip())
    save()
    return result


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''): h.update(block)
    return h.hexdigest()


def nonempty(value):
    return value not in ('', '-', 'none', 'NA')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('/work/cyu/annotation/scotti-annotation'))
    p.add_argument('--candidate', type=Path)
    p.add_argument('--out', type=Path)
    p.add_argument('--expected-proteins', type=int, default=11464)
    p.add_argument('--expected-genes', type=int, default=9234)
    a = p.parse_args()
    root = a.root
    candidate = a.candidate or root/'results/filtered_candidates/Scotti_filtered7_v1'
    out = a.out or candidate/'functional_evidence_v1'
    if out.exists(): raise FileExistsError(f'Output exists; inspect rather than overwrite: {out}')
    protein = candidate/'Scotti_nuclear.filtered7.proteins.fa'
    longest = candidate/'Scotti_nuclear.filtered7.longest.proteins.fa'
    gff = candidate/'Scotti_nuclear.filtered7.candidate.gff3'
    original = root/'results/braker3_etp_v1_retry2/braker.aa'
    swiss = root/'results/functional/01_swissprot'
    swiss_best = swiss/'scotti_vs_swissprot.best_hit.tsv'
    swiss_top = swiss/'scotti_vs_swissprot.top5.tsv'
    egg = root/'results/functional/02_eggnog/scotti_braker.emapper.annotations'
    ipr = root/'results/functional/03_interproscan/scotti_braker.tsv'
    sources = [protein, longest, gff, original, swiss_best, swiss_top, egg, ipr]
    for path in sources:
        if not path.is_file() or not path.stat().st_size: raise FileNotFoundError(path)
    lengths = fasta_lengths(protein); keep = set(lengths)
    primary = set(fasta_lengths(longest)); original_ids = set(fasta_lengths(original))
    if len(keep) != a.expected_proteins or len(primary) != a.expected_genes: raise ValueError('Unexpected candidate protein counts')
    if not primary <= keep or not keep <= original_ids: raise ValueError('Inconsistent protein sets')
    tx_gene = {}
    for line in gff.open():
        if line.startswith('##FASTA'): break
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) != 9 or f[2] != 'mRNA': continue
        at = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if at['ID'] in tx_gene or ',' in at['Parent']: raise ValueError('Unexpected transcript mapping')
        tx_gene[at['ID']] = at['Parent']
    if set(tx_gene) != keep or len(set(tx_gene.values())) != a.expected_genes: raise ValueError('GFF/protein ID mismatch')
    if len({tx_gene[x] for x in primary}) != a.expected_genes: raise ValueError('Primary proteins are not one per gene')
    out.mkdir(parents=True)
    swiss_data = {}
    for source, target, is_best in [
        (swiss_best, 'Scotti_filtered7.swissprot.best_hit.tsv', True),
        (swiss_top, 'Scotti_filtered7.swissprot.top5.tsv', False),
    ]:
        with source.open() as src, (out/target).open('x') as dst:
            for line in src:
                if not line.strip(): continue
                f = line.rstrip('\r\n').split('\t', 9)
                if len(f) != 10: raise ValueError(f'Unexpected Swiss-Prot columns: {source}')
                if f[0] not in original_ids: raise ValueError(f'Unknown Swiss-Prot query: {f[0]}')
                if f[0] not in keep: continue
                dst.write(line if line.endswith('\n') else line+'\n')
                if is_best:
                    if f[0] in swiss_data: raise ValueError('Duplicate best-hit query')
                    swiss_data[f[0]] = f
    egg_data = {}; header = None
    with egg.open() as src, (out/'Scotti_filtered7.emapper.annotations').open('x') as dst:
        for line in src:
            if not line.strip(): continue
            if line.startswith('#'):
                f = line.lstrip('#').rstrip('\r\n').split('\t')
                if f[0] in ('query', 'query_name'): header = f
                dst.write(line if line.endswith('\n') else line+'\n'); continue
            if header is None: raise ValueError('eggNOG header missing')
            f = line.rstrip('\r\n').split('\t')
            if len(f) != len(header): raise ValueError('Unexpected eggNOG column count')
            if f[0] not in original_ids: raise ValueError(f'Unknown eggNOG query: {f[0]}')
            if f[0] not in keep: continue
            if f[0] in egg_data: raise ValueError('Duplicate eggNOG query')
            egg_data[f[0]] = dict(zip(header, f))
            dst.write(line if line.endswith('\n') else line+'\n')
    if header is None: raise ValueError('eggNOG header missing')
    ipr_seen = set(); ipr_ids = defaultdict(set); ipr_desc = defaultdict(set)
    databases = defaultdict(set); signatures = defaultdict(set); ipr_go = defaultdict(set)
    with ipr.open() as src, (out/'Scotti_filtered7.interproscan.tsv').open('x') as dst:
        for line in src:
            if not line.strip(): continue
            if line.startswith('#'):
                dst.write(line); continue
            f = line.rstrip('\r\n').split('\t')
            if len(f) < 11: raise ValueError('Unexpected InterProScan column count')
            pid = f[0]
            if pid not in original_ids: raise ValueError(f'Unknown InterProScan query: {pid}')
            if pid not in keep: continue
            dst.write(line if line.endswith('\n') else line+'\n')
            ipr_seen.add(pid); databases[pid].add(f[3]); signatures[pid].add(f'{f[3]}:{f[4]}')
            if len(f) >= 13 and f[11].startswith('IPR'):
                ipr_ids[pid].add(f[11]); ipr_desc[pid].add(f'{f[11]}:{f[12]}')
            if len(f) >= 14 and nonempty(f[13]): ipr_go[pid].update(f[13].split('|'))
    fields = ['protein_id','gene_id','primary_isoform','length_aa',
        'swissprot_subject','swissprot_identity_pct','swissprot_query_coverage_pct',
        'swissprot_evalue','swissprot_bitscore','swissprot_title',
        'eggnog_row_present','eggnog_seed_ortholog','eggnog_evalue','eggnog_score',
        'eggnog_description','eggnog_preferred_name','eggnog_COG_category','eggnog_GOs',
        'eggnog_KEGG_ko','eggnog_PFAMs','interpro_member_databases','interpro_signatures',
        'interpro_ids','interpro_descriptions','interpro_GOs','review_status']
    def row(pid):
        s = swiss_data.get(pid); e = egg_data.get(pid, {})
        return [pid, tx_gene[pid], pid in primary, lengths[pid],
            s[1] if s else '', s[2] if s else '', s[6] if s else '',
            s[7] if s else '', s[8] if s else '', s[9] if s else '',
            pid in egg_data, e.get('seed_ortholog',''), e.get('evalue',''), e.get('score',''),
            e.get('Description',''), e.get('Preferred_name',''), e.get('COG_category',''),
            e.get('GOs',''), e.get('KEGG_ko',''), e.get('PFAMs',''),
            ';'.join(sorted(databases[pid])), ';'.join(sorted(signatures[pid])),
            ';'.join(sorted(ipr_ids[pid])), ';'.join(sorted(ipr_desc[pid])),
            '|'.join(sorted(ipr_go[pid])), 'evidence_only_not_curated_product']
    for filename, ids in [('Scotti_filtered7.all_isoforms.evidence.tsv', keep),
                          ('Scotti_filtered7.primary.evidence.tsv', primary)]:
        with (out/filename).open('x') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(fields)
            for pid in sorted(ids): w.writerow(row(pid))
    coverage = {}
    for label, ids in [('all_isoforms',keep),('primary',primary)]:
        coverage[label] = dict(proteins=len(ids), swissprot_hits=len(ids & swiss_data.keys()),
            eggnog_annotation_rows=len(ids & egg_data.keys()),
            interpro_member_matches=len(ids & ipr_seen),
            integrated_interpro_entries=sum(bool(ipr_ids[x]) for x in ids),
            any_of_three=len(ids & (swiss_data.keys() | egg_data.keys() | ipr_seen)))
    (out/'coverage_summary.json').write_text(json.dumps(coverage, indent=2)+'\n')
    (out/'README.txt').write_text(
        'Evidence integration for the filtered7 candidate; not a final annotation release.\n'
        'No searches rerun, no gene models edited, and no automatic product names assigned.\n'
        'Swiss-Prot hit and eggNOG row presence are not equivalent to validated function.\n'
        'Member-database matches may include low-complexity/non-specific signatures.\n'
        'Evidence absent from these resources is not a gene-removal criterion by itself.\n'
        'InterProScan TSV rows are filtered; its large source GFF3 is not rewritten.\n'
        'Filtered eggNOG comment metadata is copied from the original run; use coverage_summary.json for retained-query counts.\n'
        'Genome coordinates, protein IDs, and source search outputs are unchanged.\n')
    shutil.copy2(Path(__file__), out/'integrate_filtered7_evidence.py')
    with (out/'input.SHA256SUMS').open('x') as f:
        for path in sources: f.write(f'{digest(path)}  {path}\n')
    with (out/'SHA256SUMS').open('x') as f:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.name != 'SHA256SUMS': f.write(f'{digest(path)}  {path.name}\n')
    print(json.dumps(coverage, indent=2)); print('Output:',out)

if __name__ == '__main__': main()
