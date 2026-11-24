#!/usr/bin/env python3
"""Inventory final FASTAs and existing BUSCO summaries. No BUSCO jobs are launched.

Read-only for all source files. Writes a new report directory only. Historical
input files are compared with selected inputs by normalized ID/sequence hashes.
A present-day match is not proof of a run-time checksum when none was recorded.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import gzip
import hashlib
import json
import re
import sys

SPECIES = ['S_solidus', 'S_cotti', 'L_intestinalis', 'D_latus',
           'Spirometra_Aus1', 'S_proliferum', 'E_multilocularis',
           'H_microstoma', 'E_granulosus', 'T_solium', 'T_multiceps']
SCORES = re.compile(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],\s*'
                    r'F:([\d.]+)%,M:([\d.]+)%,n:(\d+)')


def tsv(path, rows, fields):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def text_open(path):
    with path.open('rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if compressed else path.open()


def fingerprint(path, cache):
    if not path.is_file():
        return None
    resolved = path.resolve()
    if resolved in cache:
        return cache[resolved]
    raw_hash = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            raw_hash.update(chunk)
    records = {}
    current = None
    seq_hash = None
    lengths = Counter()
    with text_open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if current is not None:
                    records[current] = seq_hash.hexdigest()
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f'Duplicate FASTA ID {current}: {path}')
                seq_hash = hashlib.sha256()
            else:
                seq = ''.join(line.split()).upper()
                if not seq:
                    continue
                if current is None:
                    raise ValueError(f'Sequence before FASTA header: {path}')
                seq_hash.update(seq.encode('ascii'))
                lengths[current] += len(seq)
    if current is not None:
        if current in records:
            raise ValueError(f'Duplicate FASTA ID {current}: {path}')
        records[current] = seq_hash.hexdigest()
    if not records or any(lengths[r] == 0 for r in records):
        raise ValueError(f'Empty FASTA or empty sequence: {path}')
    canonical = hashlib.sha256()
    for name, digest in sorted(records.items()):
        canonical.update(name.encode() + b'\t' + digest.encode() + b'\n')
    result = dict(raw_sha256=raw_hash.hexdigest(), id_sequence_sha256=canonical.hexdigest(),
                  records=len(records), residues=sum(lengths.values()))
    cache[resolved] = result
    return result


def parse_summary(path):
    text = path.read_text(errors='replace')
    match = SCORES.search(text)
    if not match:
        return None
    def capture(pattern):
        hit = re.search(pattern, text, re.M | re.I)
        return hit.group(1).strip() if hit else ''
    mode_raw = capture(r'^#?\s*BUSCO was run in mode:\s*(.+)$')
    mode = ('proteins' if mode_raw.startswith('protein') else
            'genome' if 'genome' in mode_raw else mode_raw)
    input_name = capture(r'^#?\s*Summarized benchmarking.*?for file\s+(.+)$')
    version = capture(r'^#?\s*BUSCO version is:\s*(.+)$')
    lineage = capture(r'^#?\s*The lineage dataset is:\s*(\S+)')
    creation = capture(r'Creation date:\s*([^,)]+)')
    row = dict(summary=str(path), mode=mode, pipeline=mode_raw,
               busco_version=version, lineage=lineage, dataset_date=creation,
               reported_input=input_name, **dict(zip(['C','S','D','F','M','n'], match.groups())))
    # Recover relative input paths from the run parameters only when available.
    for json_path in [path.with_suffix('.json')]:
        if json_path.is_file():
            try:
                info = json.loads(json_path.read_text())
                params = info.get('parameters', {})
                if not row['reported_input']:
                    row['reported_input'] = str(params.get('in', ''))
            except (ValueError, AttributeError):
                pass
    return row


def infer_species(row):
    haystack = row['summary'] + ' ' + row['reported_input']
    found = [sp for sp in SPECIES if re.search(
        r'(?<![A-Za-z0-9])' + re.escape(sp) + r'(?![A-Za-z0-9])', haystack)]
    if len(found) == 1:
        return found[0]
    lowered = haystack.lower()
    if '/scotti-annotation/' in lowered or '/worm/scotti/' in lowered or '/scotti_release_' in lowered:
        return 'S_cotti'
    if '/ssol-annotation/' in lowered or '/ssolidus_release_' in lowered or 'ssolidus_nuclear' in lowered:
        return 'S_solidus'
    return ''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/work/cyu/Schistocephalus_comparative_genomics'))
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, help='Optional explicit TSV: species, mode, selected_input')
    parser.add_argument('--extra-root', action='append', type=Path, default=[])
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    std = args.root/'02_standardized_data'
    ssolidus = Path('/work/shared/Ssolidus_release_v1.0_preNCBI_20260831')
    scotti = Path('/work/shared/Scotti_release_v1.0_preNCBI_20261118')
    selected = []
    for sp in SPECIES:
        for mode in ['genome', 'proteins']:
            file = std/('genomes' if mode == 'genome' else 'proteins_primary')/(sp + ('.fasta' if mode == 'genome' else '.faa'))
            if sp in ['S_solidus', 'S_cotti']:
                release, name = (ssolidus, 'Ssolidus') if sp == 'S_solidus' else (scotti, 'Scotti')
                file = release/('assembly/'+name+'_nuclear.fa' if mode == 'genome' else
                               'proteins/'+name+'_proteins.longest_isoform.v1.0.fa')
            selected.append(dict(species=sp, mode=mode, selected_input=str(file)))
    if args.inputs:
        selected = list(csv.DictReader(args.inputs.open(), delimiter='\t'))
        expected = {(sp, mode) for sp in SPECIES for mode in ['genome', 'proteins']}
        if len(selected) != 22 or {(r['species'], r['mode']) for r in selected} != expected:
            raise ValueError('Explicit input TSV must contain exactly 22 unique species/mode rows.')
    cache = {}
    for row in selected:
        print('INPUT:', row['species'], row['mode'], flush=True)
        try:
            result = fingerprint(Path(row['selected_input']), cache)
            row.update(result or {})
            row['input_status'] = 'present' if result else 'MISSING'
        except (OSError, ValueError) as error:
            row['input_status'] = 'ERROR: ' + str(error)
    tsv(args.out/'selected_inputs.tsv', selected,
        ['species','mode','selected_input','input_status','records','residues','raw_sha256','id_sequence_sha256'])
    roots = [args.root/'03_qc', ssolidus, scotti,
             Path('/work/cyu/annotation/scotti-annotation/results/annotation_qc'),
             Path('/work/cyu/annotation/scotti-annotation/results/filtered_candidates/Scotti_filtered7_v1'),
             Path('/work/cyu/annotation/ssol-annotation/results'),
             Path('/work/cyu/assembly/worm/scotti/qc'),
             Path('/work/cyu/assembly/worm/qc'), *args.extra_root]
    candidates = set()
    for root in roots:
        if root.is_dir():
            candidates.update(root.rglob('short_summary*.txt'))
    rows = []
    lookup = {(r['species'], r['mode']): r for r in selected}
    for path in sorted(candidates):
        try:
            row = parse_summary(path)
            if not row:
                continue
            row['species'] = infer_species(row)
            target = lookup.get((row['species'], row['mode']))
            row['match_status'] = 'species_or_mode_unresolved'
            if target:
                row['selected_input'] = target['selected_input']
                reported = Path(row['reported_input']) if row['reported_input'] else None
                current = fingerprint(reported, cache) if reported and reported.is_absolute() else None
                if target['input_status'] != 'present':
                    row['match_status'] = 'selected_input_missing_or_invalid'
                elif not current:
                    row['match_status'] = 'reported_input_unavailable'
                elif current['id_sequence_sha256'] == target['id_sequence_sha256']:
                    row['match_status'] = 'current_ID_sequences_match_run_time_hash_not_proven'
                else:
                    row['match_status'] = 'DIFFERENT_INPUT'
                row['current_reported_input_sha256'] = current['raw_sha256'] if current else ''
            rows.append(row)
        except (OSError, ValueError) as error:
            rows.append(dict(summary=str(path), match_status='ERROR: '+str(error)))
    fields = ['species','mode','pipeline','lineage','dataset_date','busco_version',
              'C','S','D','F','M','n','match_status','reported_input',
              'selected_input','current_reported_input_sha256','summary']
    tsv(args.out/'all_existing_busco_reports.tsv', rows, fields)
    review = []
    for target in selected:
        matches = [r for r in rows if r.get('species') == target['species']
                   and r.get('mode') == target['mode'] and r.get('lineage') == 'metazoa_odb12'
                   and r.get('match_status') == 'current_ID_sequences_match_run_time_hash_not_proven']
        signatures = {(r['C'],r['S'],r['D'],r['F'],r['M'],r['n'],r['dataset_date'],r['busco_version']) for r in matches}
        state = ('NO_MATCHING_ODB12_REPORT' if not matches else
                 'CONFLICTING_MATCHING_REPORTS' if len(signatures)>1 else
                 'matched_current_input_check_run_time_provenance')
        record = {**target, 'review_status':state, 'matching_reports':len(matches)}
        # Do not choose a report by newest date or first filename when results conflict.
        if len(signatures)==1:
            record.update({key:matches[0][key] for key in ['C','S','D','F','M','n','dataset_date','busco_version']})
        review.append(record)
    tsv(args.out/'eleven_species.genome_and_protein.review.tsv', review,
        ['species','mode','input_status','review_status','matching_reports',
         'C','S','D','F','M','n','dataset_date','busco_version','selected_input'])
    (args.out/'README.txt').write_text(
        'These are audit outputs, not a newly completed BUSCO comparison.\n'
        'Protein inputs use one representative per gene; focal species use the selected shared releases.\n'
        'Public background genomes are the standardized inputs; nuclear-only scope is not inferred.\n'
        'A current input match does not reconstruct a missing run-time hash. Check frozen provenance.\n'
        'Do not use raw primary, all-isoform, without32, odb10, or different-input runs as final odb12 results.\n'
        'No new BUSCO jobs were launched and no source files were modified.\n')
    print('\n=== 22 selected species/mode assessments ===')
    for row in review:
        print(row['species'], row['mode'], row['review_status'], sep='\t')
    print('\nReports:', args.out)


if __name__ == '__main__':
    main()
