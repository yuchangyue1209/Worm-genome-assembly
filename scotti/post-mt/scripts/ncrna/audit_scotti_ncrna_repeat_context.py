#!/usr/bin/env python3
"""Audit RepeatMasker context of the combined ncRNA candidates without removing loci.

Requires build_scotti_ncrna_candidate.py in the same directory (standard library only).
Small-RNA matches are reported separately from TE-class annotations.
"""
import argparse
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from urllib.parse import unquote

import build_scotti_ncrna_candidate as common


CATEGORIES = ('TE_class', 'small_RNA', 'unclassified', 'satellite',
              'simple_repeat', 'low_complexity', 'other')


def category(label):
    base = label.split('/')[0].rstrip('?')
    if base in {'SINE', 'LINE', 'LTR', 'DNA', 'RC', 'Retroposon', 'PLE'}:
        return 'TE_class'
    if base in {'RNA', 'rRNA', 'tRNA', 'snRNA', 'scRNA', 'srpRNA', 'snoRNA', 'Small_RNA'}:
        return 'small_RNA'
    if base in {'Unknown', 'Unclassified'}:
        return 'unclassified'
    return {'Satellite': 'satellite', 'Simple_repeat': 'simple_repeat',
            'Low_complexity': 'low_complexity'}.get(base, 'other')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('/work/cyu/annotation/scotti-annotation'))
    p.add_argument('--genome', type=Path, default=Path('/work/cyu/assembly/worm/scotti/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa'))
    p.add_argument('--repeatmasker', type=Path, default=Path('/work/cyu/assembly/worm/scotti/repeats/02_repeatmasker/Scotti_nuclear.mt_filtered.ragtag.provisional.fa.out'))
    p.add_argument('--out', type=Path)
    a = p.parse_args()
    base = a.root/'results/ncrna_v1'
    source = base/'05_ncrna_candidate_v1'
    candidate = source/'Scotti_nuclear.ncRNA.combined.candidate.gff3'
    summary_file = source/'summary.json'
    out = a.out or base/'06_ncrna_repeat_context_v1'
    common.require(not out.exists(), 'Output exists; choose a new --out: '+str(out))
    for path in [a.genome, a.repeatmasker, candidate, summary_file]:
        common.require(path.is_file() and path.stat().st_size > 0, 'Missing input: '+str(path))
    summary_in = json.loads(summary_file.read_text())
    sizes = common.fasta_lengths(a.genome)
    features = common.read_gff(candidate, sizes)
    common.require(len(features) == summary_in['combined_candidate_features'], 'Candidate count mismatch')
    wanted = {f[0] for f in features}
    intervals = defaultdict(list)
    parsed = 0
    labels = Counter()
    relevant = 0
    with a.repeatmasker.open() as handle:
        for lineno, line in enumerate(handle, 1):
            f = line.split()
            if not f or not f[0].isdigit():
                continue
            common.require(len(f) >= 15, f'Malformed RepeatMasker record at line {lineno}')
            parsed += 1
            seq, start, end = f[4], int(f[5]), int(f[6])
            common.require(seq in sizes and 1 <= start <= end <= sizes[seq],
                           f'RepeatMasker/genome coordinate mismatch at line {lineno}')
            group = category(f[10])
            labels[(f[10], group)] += 1
            if seq not in wanted:
                continue
            relevant += 1
            intervals[(seq, group)].append((start-1, end))
            intervals[(seq, 'all')].append((start-1, end))
    common.require(parsed > 0, 'No RepeatMasker records parsed')
    indexes = {key: common.Coverage(iv) for key, iv in intervals.items()}
    rows = []
    ids = set()
    for f in features:
        attrs = {k: unquote(v) for item in f[8].split(';') if '=' in item
                 for k, v in [item.split('=', 1)]}
        fid = attrs.get('ID', '')
        common.require(fid and fid not in ids, 'Missing or duplicate candidate ID')
        ids.add(fid)
        start, end = int(f[3])-1, int(f[4])
        span = end-start
        family = 'tRNA_high_confidence' if f[2] == 'tRNA' else attrs.get('Name', 'unknown')
        row = dict(feature_id=fid, family=family, feature_type=f[2], seqid=f[0],
                   start=start+1, end=end, strand=f[6], genomic_span_bp=span)
        for group in ('all',)+CATEGORIES:
            index = indexes.get((f[0], group))
            bp = index.overlap(start, end) if index else 0
            common.require(0 <= bp <= span, 'Invalid overlap length')
            row[group+'_union_bp'] = bp
            row[group+'_union_pct'] = round(100*bp/span, 3)
        # Descriptive review flags only; no exclusion or functional inference.
        row['TE_class_ge50pct'] = row['TE_class_union_bp']*2 >= span
        row['TE_class_ge80pct'] = row['TE_class_union_bp']*5 >= span*4
        row['small_RNA_ge80pct'] = row['small_RNA_union_bp']*5 >= span*4
        row['review_note'] = ('high_TE_class_overlap_review_only' if row['TE_class_ge80pct']
                              else 'no_high_TE_class_overlap_flag')
        rows.append(row)
    families = []
    for family in sorted({r['family'] for r in rows}):
        subset = [r for r in rows if r['family'] == family]
        families.append(dict(family=family, candidates=len(subset),
                             TE_class_ge50pct=sum(r['TE_class_ge50pct'] for r in subset),
                             TE_class_ge80pct=sum(r['TE_class_ge80pct'] for r in subset),
                             small_RNA_ge80pct=sum(r['small_RNA_ge80pct'] for r in subset),
                             median_TE_class_pct=round(median(r['TE_class_union_pct'] for r in subset), 3),
                             median_small_RNA_pct=round(median(r['small_RNA_union_pct'] for r in subset), 3)))
    summary = dict(status='repeat_context_audit_only_no_candidates_removed',
                   candidate_features=len(rows),
                   tRNA_features=sum(r['feature_type'] == 'tRNA' for r in rows),
                   Rfam_features=sum(r['feature_type'] != 'tRNA' for r in rows),
                   RepeatMasker_records=parsed,
                   TE_class_ge50pct=sum(r['TE_class_ge50pct'] for r in rows),
                   TE_class_ge80pct=sum(r['TE_class_ge80pct'] for r in rows),
                   small_RNA_ge80pct=sum(r['small_RNA_ge80pct'] for r in rows),
                   U1=next((r for r in families if r['family'] == 'U1'), None),
                   candidate_features_removed=0, protein_models_removed=0)
    out.mkdir(parents=True)
    common.write_tsv(out/'ncRNA_repeat_context.by_feature.tsv', list(rows[0]), rows)
    common.write_tsv(out/'ncRNA_repeat_context.by_family.tsv', list(families[0]), families)
    common.write_tsv(out/'RepeatMasker_class_mapping.tsv', ['RepeatMasker_class', 'audit_category', 'records'],
                     [dict(RepeatMasker_class=label, audit_category=group, records=n)
                      for (label, group), n in sorted(labels.items())])
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (out/'README.txt').write_text(
        'Context audit of the existing combined ncRNA candidates. No features, sequences or protein models are removed.\n'
        'RepeatMasker coordinates are 1-based closed; overlaps are calculated using merged 0-based half-open intervals.\n'
        'Both strands of RepeatMasker annotations contribute to genomic repeat context.\n'
        'Each category and all-category overlap use union-covered bases; categories can overlap each other.\n'
        'TE_class means a RepeatMasker library class annotation (SINE/LINE/LTR/DNA/RC/Retroposon/PLE), not proven TE activity.\n'
        'Small-RNA classifications, Unknown, satellite, simple repeat and low complexity are separate from TE_class.\n'
        'Unrecognized classes are retained as other and listed explicitly in RepeatMasker_class_mapping.tsv.\n'
        '50%/80% flags are descriptive project review thresholds, not validated pseudogene or deletion rules.\n'
        'Small-RNA masking alone cannot distinguish functional ncRNAs from pseudogenes. No-overlap does not prove function.\n'
        'Fractions use the genomic feature span, including introns for precursor tRNAs; mature RNA exons are not reconstructed.\n'
        'Reference: https://www.repeatmasker.org/webrepeatmaskerhelp.html\n')
    sources = [Path(__file__), Path(common.__file__)]
    for path in sources:
        shutil.copy2(path, out/path.name)
    with (out/'input.SHA256SUMS').open('x') as handle:
        for path in [candidate, summary_file, a.repeatmasker, a.genome]+sources:
            handle.write(f'{common.sha(path)}  {path}\n')
    with (out/'SHA256SUMS').open('x') as handle:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.name != 'SHA256SUMS':
                handle.write(f'{common.sha(path)}  {path.name}\n')
    # Keep the concise summary last so it is easy to copy from a terminal.
    print('Output:', out)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
