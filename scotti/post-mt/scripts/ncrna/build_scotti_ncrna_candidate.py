#!/usr/bin/env python3
"""Build a reviewable nuclear ncRNA candidate package; never edit source models.

Python standard library only. Rfam candidates are sequence-homology predictions,
not confirmed functional RNA genes. All raw Rfam hits remain in the audit table.
"""
import argparse
import csv
import hashlib
import json
import math
import shutil
from bisect import bisect_right
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path
from urllib.parse import quote, unquote


# Explicit project policy, not a universal Rfam taxonomic classifier.
POLICY = {
    'RF00001': ('5S_rRNA', 'rRNA'),
    'RF00002': ('5_8S_rRNA', 'rRNA'),
    'RF01960': ('SSU_rRNA_eukarya', 'rRNA'),
    'RF02543': ('LSU_rRNA_eukarya', 'rRNA'),
    'RF00003': ('U1', 'snRNA'), 'RF00004': ('U2', 'snRNA'),
    'RF00015': ('U4', 'snRNA'), 'RF00020': ('U5', 'snRNA'),
    'RF00026': ('U6', 'snRNA'), 'RF00548': ('U11', 'snRNA'),
    'RF00007': ('U12', 'snRNA'), 'RF00618': ('U4atac', 'snRNA'),
    'RF00619': ('U6atac', 'snRNA'),
    'RF00012': ('U3', 'snoRNA'), 'RF00049': ('SNORD36', 'snoRNA'),
    'RF00574': ('SNORD69', 'snoRNA'), 'RF00600': ('SNORA79', 'snoRNA'),
    'RF00157': ('SNORD39', 'snoRNA'),
    'RF00017': ('Metazoa_SRP', 'ncRNA'),
    'RF00009': ('RNaseP_nuc', 'ncRNA'),
    'RF00832': ('mir-71', 'ncRNA'),
}
TRNA = {'RF00005', 'RF01852'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta_lengths(path):
    result = {}
    name = None
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                require(name not in result, 'Duplicate FASTA ID: ' + name)
                result[name] = 0
            elif line.strip():
                require(name is not None, 'Sequence before FASTA header')
                result[name] += len(''.join(line.split()))
    require(bool(result) and all(result.values()), 'Empty FASTA sequence')
    return result


def read_gff(path, sizes):
    rows = []
    with path.open() as handle:
        require(handle.readline().strip() == '##gff-version 3',
                'Missing GFF3 version header: ' + str(path))
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\r\n').split('\t')
            require(len(f) == 9, 'Malformed GFF3: ' + str(path))
            require(f[0] in sizes and 1 <= int(f[3]) <= int(f[4]) <= sizes[f[0]],
                    'GFF3/genome coordinate mismatch')
            rows.append(f)
    return rows


class Coverage:
    def __init__(self, intervals):
        merged = []
        for s, e in sorted(intervals):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
            else:
                merged.append((s, e))
        self.starts = [s for s, e in merged]
        self.ends = [e for s, e in merged]
        self.prefix = [0]
        for s, e in merged:
            self.prefix.append(self.prefix[-1] + e - s)

    def before(self, x):
        i = bisect_right(self.ends, x)
        n = self.prefix[i]
        if i < len(self.starts) and self.starts[i] < x:
            n += x - self.starts[i]
        return n

    def overlap(self, s, e):
        return self.before(e) - self.before(s)


def encode(values):
    return ';'.join(k + '=' + quote(str(v), safe='._:-') for k, v in values.items())


def write_tsv(path, fields, rows):
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_gff(path, rows, sizes):
    ids = set()
    for f in rows:
        attributes = dict(item.split('=', 1) for item in f[8].split(';') if '=' in item)
        fid = unquote(attributes.get('ID', ''))
        require(fid and fid not in ids, 'Missing or duplicate output feature ID: ' + fid)
        ids.add(fid)
    with path.open('x') as handle:
        handle.write('##gff-version 3\n')
        for seq in sorted({f[0] for f in rows}):
            handle.write(f'##sequence-region {seq} 1 {sizes[seq]}\n')
        for f in sorted(rows, key=lambda f: (f[0], int(f[3]), int(f[4]), f[6], f[8])):
            handle.write('\t'.join(f) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('/work/cyu/annotation/scotti-annotation'))
    p.add_argument('--genome', type=Path, default=Path('/work/cyu/assembly/worm/scotti/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa'))
    p.add_argument('--out', type=Path)
    p.add_argument('--min-model-coverage', type=float, default=0.80)
    p.add_argument('--expected-trna', type=int, default=344)
    a = p.parse_args()
    require(0 < a.min_model_coverage <= 1, 'Invalid model coverage threshold')
    base = a.root / 'results/ncrna_v1'
    raw = base / '02_rfam/Scotti_nuclear.rfam.tblout'
    trna = base / '04_trna_candidate_v1/Scotti_nuclear.tRNA.high_confidence.candidate.gff3'
    protein = a.root / 'results/filtered_candidates/Scotti_filtered7_v1/Scotti_nuclear.filtered7.candidate.gff3'
    out = a.out or base / '05_ncrna_candidate_v1'
    require(not out.exists(), 'Output exists; preserve it and choose a new --out: ' + str(out))
    for path in [a.genome, raw, trna, protein]:
        require(path.is_file() and path.stat().st_size > 0, 'Missing input: ' + str(path))
    sizes = fasta_lengths(a.genome)
    require('Scotti_mt' not in sizes, 'Expected nuclear-only genome')
    trna_rows = read_gff(trna, sizes)
    require(len(trna_rows) == a.expected_trna and all(f[2] == 'tRNA' for f in trna_rows),
            'Unexpected high-confidence tRNA feature count or type')
    intervals = defaultdict(list)
    for f in read_gff(protein, sizes):
        if f[2] == 'CDS':
            intervals[(f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    require(bool(intervals), 'No CDS features found in protein annotation')
    cds = {k: Coverage(v) for k, v in intervals.items()}
    ti = defaultdict(list)
    for f in trna_rows:
        ti[(f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    trna_index = {k: Coverage(v) for k, v in ti.items()}
    hits = []
    local_ids = set()
    with raw.open() as handle:
        for lineno, line in enumerate(handle, 1):
            if line.startswith('#') or not line.strip():
                continue
            f = line.split(maxsplit=28)
            require(len(f) >= 28 and f[0].isdigit() and f[2].startswith('RF'),
                    f'Expected cmscan fmt 2 with model/query lengths at line {lineno}')
            seq, strand = f[3], f[11]
            start, end = sorted([int(f[9]), int(f[10])])
            ml, ql = int(f[26]), int(f[27])
            require(seq in sizes and ql == sizes[seq] and 1 <= start <= end <= ql,
                    f'Rfam/genome mismatch at line {lineno}')
            require(strand in {'+', '-'} and f[19] in {'*', '^', '$', '='}, 'Unexpected strand/overlap mark')
            require((strand == '+' and int(f[9]) <= int(f[10])) or
                    (strand == '-' and int(f[9]) >= int(f[10])), 'Strand/coordinate mismatch')
            mf, mt = sorted([int(f[7]), int(f[8])])
            require(1 <= mf <= mt <= ml, 'Invalid model coordinates')
            key = (seq, f[0])
            require(key not in local_ids, 'Duplicate query/hit index')
            local_ids.add(key)
            score, ev = float(f[16]), Decimal(f[17])
            require(math.isfinite(score) and ev.is_finite() and ev >= 0, 'Invalid score/E-value')
            hit = dict(hit_id=f'Scotti_Rfam_hit_{len(hits)+1:06d}', query_hit_index=f[0],
                       seqid=seq, start=start, end=end, strand=strand, family=f[1], accession=f[2],
                       clan=f[5], model_from=mf, model_to=mt, model_length=ml,
                       model_coverage=(mt-mf+1)/ml, trunc=f[12], score=score,
                       evalue=f[17], inclusion=f[18], original_olp=f[19],
                       better_overlap_query_index=f[20], original_winner_query_index=f[23],
                       description=f[28] if len(f) > 28 else '', feature_id='',
                       feature_type=POLICY.get(f[2], ('', ''))[1], decision='', reasons='')
            require(f[2] not in POLICY or f[1] == POLICY[f[2]][0],
                    'Accession/name policy mismatch: ' + f[2] + ' ' + f[1])
            index = cds.get((seq, strand))
            hit['same_strand_CDS_overlap_bp'] = index.overlap(start-1, end) if index else 0
            index = trna_index.get((seq, strand))
            hit['same_strand_high_confidence_tRNA_overlap_bp'] = index.overlap(start-1, end) if index else 0
            reasons = []
            if f[2] in TRNA:
                hit['decision'] = 'crosscheck_only_tRNAscan_is_primary'
                reasons.append('Rfam_tRNA_not_added_to_344_tRNAscan_candidates')
            else:
                if f[18] != '!': reasons.append('not_included_by_original_scan')
                if f[19] == '=': reasons.append('lower_scoring_overlap_Rfam_documented_filter')
                if f[2] not in POLICY: reasons.append('family_outside_project_core_policy')
                if f[12] != 'no': reasons.append('truncated_hit')
                if hit['model_coverage'] < a.min_model_coverage: reasons.append('short_model_span')
                if hit['same_strand_high_confidence_tRNA_overlap_bp']: reasons.append('overlaps_curated_tRNA_same_strand')
                hit['decision'] = 'review_only' if reasons else 'core_candidate'
            hit['reasons'] = ';'.join(reasons)
            hits.append(hit)
    require(bool(hits), 'No Rfam records')

    # Resolve residual overlaps among eligible candidates, including tied hits.
    # Opposite-strand hits remain separate. All displaced hits stay in the audit.
    selected = defaultdict(list)
    for h in sorted((h for h in hits if h['decision'] == 'core_candidate'),
                    key=lambda h: (Decimal(h['evalue']), -h['score'], h['hit_id'])):
        key = (h['seqid'], h['strand'])
        conflicts = [old for old in selected[key]
                     if old['start'] <= h['end'] and h['start'] <= old['end']]
        if conflicts:
            h['decision'] = 'review_only'
            h['reasons'] = 'residual_same_strand_candidate_overlap:' + ','.join(old['hit_id'] for old in conflicts)
        else:
            selected[key].append(h)
    core = sorted((h for h in hits if h['decision'] == 'core_candidate'),
                  key=lambda h: (h['seqid'], h['start'], h['end'], h['strand'], h['accession']))
    rfam_rows = []
    for i, h in enumerate(core, 1):
        h['feature_id'] = f'Scotti_Rfam_{i:05d}'
        note = 'sequence_homology_candidate_function_and_pseudogene_status_not_validated'
        if h['accession'] == 'RF00832':
            note += ';miRNA_precursor_homology_only_mature_miRNA_not_defined'
        attr = encode(dict(ID=h['feature_id'], Name=h['family'], Dbxref='Rfam:'+h['accession'],
                           model_coverage=f"{h['model_coverage']:.4f}", source_hit=h['hit_id'],
                           original_overlap_mark=h['original_olp'], Note=note))
        rfam_rows.append([h['seqid'], 'Infernal', h['feature_type'], str(h['start']),
                          str(h['end']), str(h['score']), h['strand'], '.', attr])
    combined = trna_rows + rfam_rows
    summary = dict(status='nuclear_ncRNA_candidate_not_final_functional_gene_catalogue',
                   high_confidence_tRNA_features=len(trna_rows), raw_Rfam_hits=len(hits),
                   Rfam_core_candidate_features=len(core), combined_candidate_features=len(combined),
                   Rfam_decisions=dict(Counter(h['decision'] for h in hits)),
                   Rfam_core_families=dict(sorted(Counter(h['family'] for h in core).items())),
                   Rfam_core_candidates_overlapping_same_strand_CDS=sum(h['same_strand_CDS_overlap_bp'] > 0 for h in core),
                   min_model_span_fraction=a.min_model_coverage,
                   protein_models_removed=0,
                   notes=['80% model-span default is an explicit project heuristic, not an Rfam requirement.',
                          'Model span is not alignment percent identity or proof of a complete functional transcript.',
                          'Rfam E-values are preserved as reported; no genome-wide rescaling was performed.',
                          'Non-core family models and overlapping/partial hits remain available for review.'])
    out.mkdir(parents=True)
    write_gff(out/'Scotti_nuclear.Rfam.core.candidate.gff3', rfam_rows, sizes)
    write_gff(out/'Scotti_nuclear.ncRNA.combined.candidate.gff3', combined, sizes)
    fields = list(hits[0])
    write_tsv(out/'Rfam_all_hits.review.tsv', fields, hits)
    write_tsv(out/'Rfam_core_candidates.tsv', fields, core)
    family_rows = []
    for family in sorted({h['family'] for h in hits}):
        subset = [h for h in hits if h['family'] == family]
        counts = Counter(h['decision'] for h in subset)
        family_rows.append(dict(family=family, accession=subset[0]['accession'], raw_hits=len(subset),
                                core_candidates=counts['core_candidate'], review_only=counts['review_only'],
                                tRNA_crosscheck_only=counts['crosscheck_only_tRNAscan_is_primary']))
    write_tsv(out/'Rfam_selection_summary.tsv', list(family_rows[0]), family_rows)
    write_tsv(out/'core_family_policy.tsv', ['accession', 'family', 'GFF_type'],
              [dict(accession=acc, family=name, GFF_type=kind) for acc, (name, kind) in sorted(POLICY.items())])
    with (out/'Scotti_nuclear.ncRNA.combined.candidate.bed').open('x') as handle:
        for f in sorted(combined, key=lambda f: (f[0], int(f[3]), int(f[4]))):
            fid = next(v[3:] for v in f[8].split(';') if v.startswith('ID='))
            handle.write('\t'.join([f[0], str(int(f[3])-1), f[4], fid, '.', f[6]])+'\n')
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (out/'README.txt').write_text(
        'Candidate-only nuclear ncRNA package. Protein GFF3, FASTA and all source scans are unchanged.\n'
        'tRNA: preserve the 344 corrected tRNAscan-SE high-confidence precursor loci, including the predicted SeC.\n'
        'Rfam: included original hits, explicit project family policy, nontruncated, model span >= configured threshold,\n'
        'excluding original overlap mark =, excluding same-strand overlap with curated tRNA. Marks * ^ $ are eligible.\n'
        'Residual same-strand candidate overlaps: lowest reported E-value, then highest bit score, then stable hit ID.\n'
        'This is a conservative subset, not an exhaustive catalogue or a claim that rejected hits are false positives.\n'
        'All raw hits and processing reasons are retained. Bacterial/archaeal model names alone do not prove contamination.\n'
        'Candidate RNA homology can include pseudogenes/repeat-derived copies; high copy counts are not functional gene counts.\n'
        'mir-71 is exported as generic ncRNA precursor homology; no mature-miRNA coordinates or expression validation.\n'
        'Flat features only. GFF3 uses 1-based closed coordinates; BED uses 0-based half-open coordinates.\n'
        'CDS overlap is an audit, not a reason to remove protein models. External GFF3 validation must be run separately.\n'
        'Reference: https://docs.rfam.org/genome-annotation.html\n')
    shutil.copy2(Path(__file__), out/Path(__file__).name)
    with (out/'input.SHA256SUMS').open('x') as handle:
        for path in [a.genome, raw, trna, protein]:
            handle.write(f'{sha(path)}  {path}\n')
    with (out/'SHA256SUMS').open('x') as handle:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.name != 'SHA256SUMS':
                handle.write(f'{sha(path)}  {path.name}\n')
    print(json.dumps(summary, indent=2))
    print('Output:', out)


if __name__ == '__main__':
    main()
