#!/usr/bin/env python3
"""Fill missing genome-mode BUSCO assessments from the reviewed 22-input audit.

Requires audit_final_busco_inputs.py in the same directory. Reuses all protein
results and the matched solidus genome result. Sources are never overwritten.
"""
from pathlib import Path
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys

from audit_final_busco_inputs import SPECIES, fingerprint, parse_summary, tsv

MATCHED = 'matched_current_input_check_run_time_provenance'
DATE = '2026-05-22'
VERSION = '6.0.0'


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for data in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(data)
    return digest.hexdigest()


def database_files(folder):
    files = [p for p in folder.iterdir() if p.is_file() and
             (p.name in ['dataset.cfg','scores_cutoff','lengths_cutoff'] or p.name.startswith('refseq'))]
    files.extend(p for p in (folder/'hmms').glob('*') if p.is_file())
    return {str(p.relative_to(folder)):sha(p) for p in sorted(files)}


def dataset_problems(folder):
    problems = []
    if folder.name != 'metazoa_odb12':
        problems.append('directory name must be metazoa_odb12')
    cfg = folder/'dataset.cfg'
    if not cfg.is_file():
        return problems + ['missing dataset.cfg']
    try:
        text = cfg.read_text()
    except OSError as error:
        return problems + ['cannot read dataset.cfg: '+str(error)]
    if not re.search(r'(?m)^\s*creation_date\s*=\s*'+re.escape(DATE)+r'\s*$', text):
        problems.append('creation_date must be '+DATE)
    if not re.search(r'(?mi)^\s*number_of_BUSCOs\s*=\s*672\s*$', text):
        problems.append('number_of_BUSCOs must be 672')
    if not (folder/'hmms').is_dir() or not any(p.is_file() and p.stat().st_size for p in (folder/'hmms').glob('*.hmm')):
        problems.append('missing/empty hmms/*.hmm')
    if not (folder/'scores_cutoff').is_file() or not (folder/'scores_cutoff').stat().st_size:
        problems.append('missing/empty scores_cutoff')
    if not any(p.is_file() and p.stat().st_size for p in folder.glob('refseq_db.faa*')):
        problems.append('missing/empty refseq_db.faa or refseq_db.faa.gz')
    # Original odb12 datasets intentionally have no lengths_cutoff. Do not
    # require or synthesize it. odb10/odb12.2 have different length-cutoff logic.
    # https://busco.ezlab.org/busco_userguide.html#match-thresholds
    return problems


def dataset_ok(folder):
    return not dataset_problems(folder)


def find_database(explicit, root, audit):
    if explicit:
        folder = explicit.resolve()
        problems = dataset_problems(folder)
        if problems:
            raise ValueError('Requested dataset rejected: '+str(folder)+'; '+ '; '.join(problems))
        return folder
    candidates = set()
    # Recover the actual dataset location from existing run metadata when present.
    for row in read_tsv(audit/'all_existing_busco_reports.tsv'):
        summary = Path(row.get('summary', ''))
        metadata = summary.with_suffix('.json')
        if metadata.is_file():
            try:
                info = json.loads(metadata.read_text()).get('parameters', {})
                for key in ['lineage_dataset', 'download_path']:
                    value = info.get(key)
                    if isinstance(value, str) and value:
                        candidate = Path(value).expanduser()
                        if candidate.is_absolute():
                            candidates.update([candidate, candidate/'lineages/metazoa_odb12', candidate/'metazoa_odb12'])
            except (ValueError, AttributeError):
                pass
    # Bounded project/database search. No new lineage is downloaded automatically.
    roots = [root, Path('/work/cyu/databases'), Path('/work/cyu/annotation/scotti-annotation'),
             Path('/work/cyu/annotation/ssol-annotation'), Path('/work/cyu/assembly/worm/scotti/qc'),
             Path.home()/'busco_downloads', Path.cwd()/'busco_downloads']
    for parent in roots:
        if parent.is_dir():
            candidates.update(p.parent for p in parent.rglob('dataset.cfg') if p.parent.name == 'metazoa_odb12')
    valid = sorted({p.resolve() for p in candidates if dataset_ok(p)})
    if not valid:
        for candidate in sorted(candidates):
            print('REJECTED DATASET:', candidate, '; '.join(dataset_problems(candidate)), file=sys.stderr)
        raise ValueError('No complete local metazoa_odb12 dataset dated '+DATE+' was found. Supply --lineage-db with its actual directory. No jobs launched.')
    hashes = [database_files(p) for p in valid]
    if any(row != hashes[0] for row in hashes[1:]):
        raise ValueError('Multiple local datasets have different contents. Select the intended copy using --lineage-db: '+str(valid))
    return valid[0]


def persist(path, data):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/work/cyu/Schistocephalus_comparative_genomics'))
    parser.add_argument('--audit-dir', type=Path)
    parser.add_argument('--lineage-db', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--cpu', type=int, default=20)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if args.cpu < 1:
        raise ValueError('--cpu must be positive')
    if args.resume and not args.out:
        raise ValueError('--resume requires the exact previous --out directory')
    if args.audit_dir:
        audit = args.audit_dir.resolve()
    else:
        audits = sorted(p for p in (args.root/'03_qc').glob('final_busco_audit_*')
                        if (p/'selected_inputs.tsv').is_file() and
                        (p/'eleven_species.genome_and_protein.review.tsv').is_file())
        if not audits:
            raise ValueError('No completed final_busco_audit_* directory found; provide --audit-dir')
        audit = audits[-1].resolve()
    inputs = read_tsv(audit/'selected_inputs.tsv')
    review = read_tsv(audit/'eleven_species.genome_and_protein.review.tsv')
    expected = {(sp, mode) for sp in SPECIES for mode in ['genome','proteins']}
    if len(inputs) != 22 or len(review) != 22 or {(r['species'],r['mode']) for r in inputs} != expected or {(r['species'],r['mode']) for r in review} != expected:
        raise ValueError('Expected exactly eleven species and both modes in the completed audit')
    if any(r['input_status'] != 'present' for r in inputs):
        raise ValueError('At least one selected input was missing or invalid in the audit')
    if any(r['review_status'] != MATCHED for r in review if r['mode'] == 'proteins'):
        raise ValueError('Protein results are incomplete/conflicting; resolve the audit before this genome-only batch')
    matched_genomes = [r for r in review if r['mode'] == 'genome' and r['review_status'] == MATCHED]
    if not matched_genomes:
        raise ValueError('No matched genome run is available to establish the reference pipeline')
    for row in review:
        if row['review_status'] == MATCHED and (row['busco_version'] != VERSION or row['dataset_date'] != DATE or row['n'] != '672'):
            raise ValueError('Matched reports do not share BUSCO 6.0.0 / odb12 '+DATE+' / n=672')
    reports = read_tsv(audit/'all_existing_busco_reports.tsv')
    pipelines = {r['pipeline'] for r in reports if r.get('mode') == 'genome'
                 and r.get('match_status') == 'current_ID_sequences_match_run_time_hash_not_proven'
                 and r.get('lineage') == 'metazoa_odb12'}
    if not pipelines or any('min' not in x.lower() for x in pipelines):
        raise ValueError('Existing genome pipeline is not unambiguously Miniprot; inspect pipeline column: '+str(pipelines))
    missing = [r for r in review if r['mode'] == 'genome' and r['review_status'] == 'NO_MATCHING_ODB12_REPORT']
    if any(r['mode'] == 'genome' and r['review_status'] not in [MATCHED,'NO_MATCHING_ODB12_REPORT'] for r in review):
        raise ValueError('Conflicting genome reports must be resolved before running this batch')
    busco = shutil.which('busco')
    if not busco:
        raise ValueError('busco is not available; activate genome_assembly first')
    version = subprocess.run([busco,'--version'], capture_output=True, text=True, check=True)
    if not re.search(r'(?<![\d.])6\.0\.0(?![\d.])', version.stdout+version.stderr):
        raise ValueError('Use the same BUSCO 6.0.0 installation as the matched reports')
    help_text = subprocess.run([busco,'--help'], capture_output=True, text=True, check=True).stdout
    if '--miniprot' not in help_text:
        raise ValueError('Local BUSCO help does not expose --miniprot')
    print('Checking the existing local odb12 dataset...', flush=True)
    dataset = find_database(args.lineage_db, args.root, audit)
    db_hashes = database_files(dataset)
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S')
    out = (args.out or args.root/'03_qc'/('final_genome_busco_odb12_'+timestamp)).absolute()
    run_spec = dict(audit=str(audit), busco=busco, busco_version=VERSION,
                    lineage=str(dataset), dataset_date=DATE, n=672, pipeline='miniprot',
                    database_files_SHA256=db_hashes, species_to_run=[r['species'] for r in missing])
    if args.resume:
        previous = json.loads((out/'run_spec.json').read_text())
        if previous != run_spec:
            raise ValueError('Resume metadata differ from the saved inputs/database/tool specification')
    else:
        out.mkdir(parents=True, exist_ok=False)
        (out/'inputs').mkdir()
        (out/'runs').mkdir()
        (out/'logs').mkdir()
        persist(out/'run_spec.json', run_spec)
        shutil.copyfile(audit/'selected_inputs.tsv', out/'audited_selected_inputs.tsv')
        shutil.copyfile(audit/'eleven_species.genome_and_protein.review.tsv', out/'reused_report_values.tsv')
    print('Audit:',audit, '\nPinned dataset:',dataset, '\nOutput:',out, flush=True)
    print('Genomes to assess:', ', '.join(run_spec['species_to_run']), flush=True)
    cache = {}
    for row in inputs:
        print('VERIFY INPUT:',row['species'],row['mode'],flush=True)
        observed = fingerprint(Path(row['selected_input']), cache)
        if not observed or observed['raw_sha256'] != row['raw_sha256'] or observed['id_sequence_sha256'] != row['id_sequence_sha256']:
            raise ValueError('Input changed after audit: '+row['selected_input'])
    input_lookup = {(r['species'],r['mode']):r for r in inputs}
    new_rows = {}
    for row in missing:
        sp = row['species']
        audited = input_lookup[(sp,'genome')]
        source = Path(audited['selected_input'])
        frozen = out/'inputs'/(sp+'.fa')
        if not frozen.exists():
            temporary = frozen.with_suffix('.fa.partial')
            if temporary.exists():
                raise ValueError('Incomplete staged FASTA exists; inspect before resuming: '+str(temporary))
            with source.open('rb') as handle:
                compressed = handle.read(2) == b'\x1f\x8b'
            if compressed:
                with gzip.open(source,'rb') as reader, temporary.open('xb') as writer:
                    shutil.copyfileobj(reader,writer)
            else:
                with source.open('rb') as reader, temporary.open('xb') as writer:
                    shutil.copyfileobj(reader,writer)
            if fingerprint(temporary,{})['id_sequence_sha256'] != audited['id_sequence_sha256'] or sha(source) != audited['raw_sha256']:
                raise ValueError('Source changed during staging: '+sp)
            temporary.rename(frozen)
        frozen_fp = fingerprint(frozen,{})
        if frozen_fp['id_sequence_sha256'] != audited['id_sequence_sha256']:
            raise ValueError('Frozen input differs: '+sp)
        record_file = out/'logs'/(sp+'.run.json')
        record = dict(species=sp, source_input=str(source), source_raw_sha256=audited['raw_sha256'],
                      frozen_input=str(frozen), frozen_raw_sha256=frozen_fp['raw_sha256'],
                      id_sequence_sha256=frozen_fp['id_sequence_sha256'])
        if record_file.exists() and json.loads(record_file.read_text()) != record:
            raise ValueError('Saved run/input record differs: '+sp)
        if not record_file.exists():
            persist(record_file,record)
        run_folder = out/'runs'/sp
        summary_files = list(run_folder.glob('short_summary.specific.metazoa_odb12.*.txt'))
        done = out/'logs'/(sp+'.done.json')
        if done.is_file():
            saved = json.loads(done.read_text())
            if len(summary_files)!=1 or sha(summary_files[0]) != saved['summary_sha256']:
                raise ValueError('Completed summary changed: '+sp)
            print('VERIFIED completed:',sp,flush=True)
        else:
            command = [busco,'-i',str(frozen),'-m','genome','-l',str(dataset),
                       '--offline','--miniprot','-c',str(args.cpu),
                       '--out_path',str(out/'runs'),'-o',sp]
            if run_folder.exists():
                if not args.resume:
                    raise ValueError('Unexpected existing run folder: '+sp)
                command.append('--restart')
            print('START:',sp,datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
            with (out/'logs'/(sp+'.console.log')).open('a') as log:
                log.write('\nCOMMAND: '+json.dumps(command)+'\n');log.flush()
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            summary_files = list(run_folder.glob('short_summary.specific.metazoa_odb12.*.txt'))
            if len(summary_files)!=1 or sha(frozen) != frozen_fp['raw_sha256']:
                raise ValueError('Missing/ambiguous summary or modified frozen input: '+sp)
        summary = parse_summary(summary_files[0])
        if not summary or summary['mode']!='genome' or summary['lineage']!='metazoa_odb12' or summary['n']!='672' or summary['dataset_date']!=DATE or summary['busco_version']!=VERSION:
            raise ValueError('Unexpected result metadata: '+sp)
        reported = Path(summary['reported_input'])
        if not reported.is_absolute() or reported.resolve() != frozen.resolve():
            raise ValueError('Summary is not linked to this frozen input: '+sp)
        new_rows[sp] = dict(species=sp,mode='genome',**{k:summary[k] for k in ['C','S','D','F','M','n','dataset_date','busco_version','pipeline','summary']},
                           input_provenance='frozen_input_run_hash_verified',selected_input=str(source))
        persist(done,dict(summary_sha256=sha(summary_files[0]), frozen_raw_sha256=frozen_fp['raw_sha256']))
        print('PASS:',sp,'C='+summary['C'],flush=True)
    if database_files(dataset) != db_hashes:
        raise ValueError('Pinned database core files changed during this batch')
    combined = []
    for row in review:
        if row['mode']=='genome' and row['species'] in new_rows:
            combined.append(new_rows[row['species']])
        else:
            combined.append({**row,'input_provenance':'existing_report_current_input_match_run_time_hash_not_recorded'})
    fields=['species','mode','C','S','D','F','M','n','dataset_date','busco_version','pipeline','input_provenance','selected_input','summary']
    combined_file=out/'eleven_species.genome_and_protein.odb12.tsv'
    if combined_file.exists():
        combined_file.unlink()  # Only this script's reproducible summary inside its own run directory.
    tsv(combined_file,combined,fields)
    persist(out/'COMPLETE.json',dict(status='completed_22_assessments',new_genome_runs=len(missing),
                                    reused_existing_assessments=22-len(missing),comparison=str(combined_file)))
    print('COMPLETE:',combined_file,flush=True)


if __name__=='__main__':
    try:
        main()
    except Exception as error:
        print('STOPPED:',error, file=sys.stderr)
        print('Existing source data and completed outputs are preserved.', file=sys.stderr)
        sys.exit(1)
