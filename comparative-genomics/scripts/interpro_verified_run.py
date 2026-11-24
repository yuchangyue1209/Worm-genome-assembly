#!/usr/bin/env python3
"""Run a fresh InterPro job on a derived FASTA; preserve existing evidence.

A nonempty TSV alone is not a completion record. Historical TSVs without our
completion manifest are never overwritten, adopted or silently rerun.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def prepare(path):
    records = {}
    identifier = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            fields = line[1:].split()
            if not fields or fields[0] in records:
                raise ValueError('Empty or duplicate FASTA ID: ' + line)
            identifier = fields[0]
            records[identifier] = ''
        elif line.strip():
            if identifier is None:
                raise ValueError('Sequence before FASTA header')
            records[identifier] += line.strip()
    if not records:
        raise ValueError('Empty FASTA')
    audit = []
    for identifier, original in records.items():
        # Only one terminal stop marker may be removed. No internal stops or
        # other substitutions are silently repaired; release FASTAs stay intact.
        sequence = original[:-1] if original.endswith('*') else original
        if not re.fullmatch('[A-Z]+', sequence):
            raise ValueError(f'Invalid/empty search sequence or internal stop: {identifier}')
        records[identifier] = sequence
        audit.append(dict(id=identifier, stored_length=len(original),
                          search_length=len(sequence), terminal_stop_removed=original.endswith('*'),
                          search_md5=hashlib.md5(sequence.encode('ascii')).hexdigest()))
    return records, audit


def validate_tsv(path, records):
    annotated = set()
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        row = line.split('\t')
        if len(row) < 11 or row[0] not in records:
            raise ValueError(f'Invalid InterPro row/ID at {path}:{number}')
        sequence = records[row[0]]
        if row[1].lower() != hashlib.md5(sequence.encode('ascii')).hexdigest() or int(row[2]) != len(sequence):
            raise ValueError(f'InterPro sequence MD5/length mismatch: {row[0]} at row {number}')
        annotated.add(row[0])
    return len(annotated)


def run(protein, prefix, executable, threads, log):
    protein, prefix, executable, log = map(Path, (protein, prefix, executable, log))
    if threads < 1:
        raise ValueError('threads must be positive')
    records, audit = prepare(protein)
    identity = dict(schema_version=1, input_path=str(protein.resolve()),
                    input_sha256=sha256(protein), executable_path=str(executable.resolve()),
                    executable_sha256=sha256(executable),
                    options=['-f', 'TSV', '-goterms', '-pa', '-cpu', str(threads)])
    result = Path(str(prefix) + '.tsv')
    search = Path(str(prefix) + '.search.faa')
    provenance = Path(str(prefix) + '.input.json')
    complete = Path(str(prefix) + '.complete.json')
    if complete.is_file():
        saved = json.loads(complete.read_text())
        if saved.get('identity') != identity or saved.get('status') != 'process_completed_and_tsv_validated':
            raise ValueError('Completion identity differs; choose a new output directory')
        for key, path in [('result_sha256', result), ('search_sha256', search), ('provenance_sha256', provenance)]:
            if not path.is_file() or sha256(path) != saved.get(key):
                raise ValueError(f'Completed output changed/missing: {path}')
        validate_tsv(result, records)
        print('VERIFIED SKIP:', result)
        return
    # This also catches interrupted jobs and stale outputs under the prefix.
    occupied = [p for p in prefix.parent.glob(prefix.name + '.*')]
    if occupied or log.exists() or log.is_symlink():
        raise FileExistsError(f'Unverified existing output; preserve it and audit or choose a new directory: {prefix}')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    with search.open('x') as handle:
        for identifier, sequence in records.items():
            handle.write(f'>{identifier}\n{sequence}\n')
    with provenance.open('x') as handle:
        json.dump(dict(identity=identity, proteins=audit), handle, indent=2)
    command = [str(executable), '-i', str(search), '-b', str(prefix), *identity['options']]
    with log.open('x') as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    if not result.is_file():
        raise FileNotFoundError(f'InterPro exited successfully but TSV is missing: {result}')
    annotated = validate_tsv(result, records)
    if sha256(protein) != identity['input_sha256'] or sha256(executable) != identity['executable_sha256']:
        raise ValueError('Input/executable changed during run; no completion record written')
    with complete.open('x') as handle:
        json.dump(dict(identity=identity, status='process_completed_and_tsv_validated',
                       result_sha256=sha256(result), search_sha256=sha256(search),
                       provenance_sha256=sha256(provenance), proteins=len(records),
                       proteins_with_tsv_rows=annotated,
                       limitation='Executable hash is recorded; database contents are not independently fingerprinted.'), handle, indent=2)
    print('PASS:', result, 'proteins with annotation rows:', annotated, '/', len(records))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protein', required=True, type=Path)
    parser.add_argument('--prefix', required=True, type=Path)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--threads', required=True, type=int)
    parser.add_argument('--log', required=True, type=Path)
    args = parser.parse_args()
    run(args.protein, args.prefix, args.executable, args.threads, args.log)
