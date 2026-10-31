#!/usr/bin/env python3
"""Paired SortMeRNA screening in bounded, sequential chunks. Standard library only.
Run on stickleback inside rrna_filter. Never edits original FASTQs.
Resume rescans the input stream but skips validated completed chunks.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime

VERSION = '1.0'
SUFFIX = re.compile(br'/[12]$')
GIB = 1024**3


def fastq(handle):
    header = handle.readline(1_000_001)
    if not header:
        return None
    seq, plus, qual = [handle.readline(1_000_001) for _ in range(3)]
    lines = (header, seq, plus, qual)
    if any(not line or len(line) > 1_000_000 for line in lines):
        raise ValueError('Truncated FASTQ record or unexpectedly long line')
    if (not header.startswith(b'@') or not plus.startswith(b'+') or
        not seq.rstrip(b'\r\n') or
        len(seq.rstrip(b'\r\n')) != len(qual.rstrip(b'\r\n'))):
        raise ValueError('Invalid FASTQ record')
    return lines


def name(record):
    return SUFFIX.sub(b'', record[0][1:].split()[0])


def paired(a, b):
    x, y = fastq(a), fastq(b)
    if x is None and y is None:
        return None
    if x is None or y is None or name(x) != name(y):
        raise ValueError('Mate IDs/counts do not match')
    return x, y


def open_reads(path):
    return gzip.open(path, 'rb') if str(path).endswith('.gz') else open(path, 'rb')


def collect_names(a, b):
    names = set()
    with open_reads(a) as x, open_reads(b) as y:
        while True:
            pair = paired(x, y)
            if pair is None:
                break
            key = name(pair[0])
            if key in names:
                raise ValueError(f'Duplicate input read ID: {key!r}')
            names.add(key)
    return names


def validate_groups(expected, groups):
    """Consumes expected set, rejecting orphan mates, duplicates and dropped IDs."""
    counts = {}
    for label, (a, b) in groups.items():
        count = 0
        with open_reads(a) as x, open_reads(b) as y:
            while True:
                pair = paired(x, y)
                if pair is None:
                    break
                key = name(pair[0])
                if key not in expected:
                    raise ValueError(f'Duplicate or unknown output ID: {key!r}')
                expected.remove(key)
                count += 1
        counts[label] = count
    if expected:
        raise ValueError(f'{len(expected)} input pairs missing from output')
    return counts


def pilot_check(qc):
    pilot = qc / '03_rrna_pilot'
    smr = qc / '04_sortmerna_pilot'
    names = collect_names(pilot/'pilot_R1.fastq', pilot/'pilot_R2.fastq')
    total = len(names)
    if total != 100000:
        raise ValueError(f'Expected 100000 pilot pairs; found {total}')
    groups = {label: (smr/f'{label}_fwd.fq', smr/f'{label}_rev.fq')
              for label in ('rrna', 'non_rrna')}
    counts = validate_groups(names, groups)
    print('PILOT PASS: paired FASTQs and input/output IDs conserved.', flush=True)
    for label, count in counts.items():
        print(f'  {label}: {count:,} pairs ({count/total:.2%})', flush=True)
    return counts


def atomic_json(path, data):
    partial = Path(str(path)+'.partial')
    with open(partial, 'w') as handle:
        json.dump(data, handle, indent=2)
        handle.write('\n')
    os.replace(partial, path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def free_space(path, minimum):
    free = shutil.disk_usage(path).free
    if free < minimum:
        raise RuntimeError(f'Insufficient free space at {path}: {free/GIB:.1f} GiB')


def output_file(prefix, mate):
    hits = [Path(str(prefix)+f'_{mate}'+ext)
            for ext in ('.fq', '.fastq', '.fq.gz', '.fastq.gz')]
    found = [p for p in hits if p.is_file()]
    if len(found) != 1:
        raise ValueError(f'Expected one {mate} output for {prefix}; found {found}')
    return found[0]


def compress(source, target):
    if target.exists():
        raise FileExistsError(f'Refusing to overwrite output: {target}')
    partial = Path(str(target)+'.partial')
    # Preserve FASTQ content; sortmerna may produce plain or gzip output.
    with open_reads(source) as src, gzip.open(partial, 'wb', compresslevel=1) as dst:
        shutil.copyfileobj(src, dst, length=4*1024*1024)
    os.replace(partial, target)


def main(a):
    binary = shutil.which('sortmerna')
    if binary is None:
        raise RuntimeError('Activate rrna_filter: sortmerna is not in PATH')
    root = a.root.expanduser().resolve()
    qc = root/'results/rnaseq_qc'
    ref = qc/'02_overrepresented_check/solidus_nuclear_rrna.fa'
    inputs = [root/'data/rnaseq/raw'/f'PCH03_802750_S1_L002_R{i}_001.fastq.gz'
              for i in (1,2)]
    for path in [ref]+inputs:
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError(path)
    pilot_check(qc)  # Always verify the successful small run before processing all reads.
    if a.pilot_only:
        return
    if a.chunk_pairs <= 0 or a.threads <= 0:
        raise ValueError('Chunk size and thread count must be positive')

    report = qc/'05_sortmerna_full'
    data = root/'data/rnaseq/rrna_screened_v1'
    scratch = a.scratch.expanduser().resolve()
    version = subprocess.run([binary,'--version'],text=True,stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,check=False).stdout
    config = dict(script_version=VERSION, chunk_pairs=a.chunk_pairs,
                  reference=str(ref), reference_sha256=sha256(ref),
                  inputs=[dict(path=str(p), size=p.stat().st_size,
                               mtime_ns=p.stat().st_mtime_ns) for p in inputs],
                  sortmerna_version=version, paired_policy='paired_in',
                  expected_input_pairs=a.expected_pairs)
    if a.resume:
        old = json.loads((report/'configuration.json').read_text())
        if old != config:
            raise ValueError('Inputs, chunk size, tool version or configuration changed; cannot resume')
        scratch.mkdir(parents=True,exist_ok=True)
    else:
        if any(p.exists() for p in (report,data,scratch)):
            raise FileExistsError('Output/work directory exists. Use --resume only for this same run.')
        report.mkdir(parents=True)
        data.mkdir(parents=True)
        scratch.mkdir(parents=True)
        atomic_json(report/'configuration.json',config)
        shutil.copy2(ref,report/'reference.solidus_nuclear_rrna.fa')
        try:
            shutil.copy2(Path(__file__),report/Path(__file__).name)
        except shutil.SameFileError:
            pass
    for label in ('rrna','non_rrna'):
        (data/label).mkdir(exist_ok=True)

    # Fixed outputs from a failed, uncheckpointed chunk are retained and require review.
    # No input FASTQ or previous biological result is deleted by this program.
    completed=[]
    total_input=0
    with gzip.open(inputs[0],'rb') as r1, gzip.open(inputs[1],'rb') as r2:
        for index in itertools.count(1):
            chunk=f'chunk{index:05d}'
            done=report/f'{chunk}.done.json'
            if done.exists():
                saved=json.loads(done.read_text())
                for item in saved['files']:
                    path=Path(item['path'])
                    if not path.is_file() or path.stat().st_size!=item['size'] or sha256(path)!=item['sha256']:
                        raise ValueError(f'Completed output changed: {path}')
                for _ in range(saved['input_pairs']):
                    if paired(r1,r2) is None:
                        raise ValueError('Input ended before checkpoint')
                total_input+=saved['input_pairs']
                completed.append(saved)
                print(f'RESUME: {chunk}, {saved["input_pairs"]:,} pairs; mapping skipped',flush=True)
                continue

            pair=paired(r1,r2)
            if pair is None:
                break
            free_space(scratch,100*GIB)
            free_space(data,100*GIB)
            work=scratch/(chunk+'.'+datetime.now().strftime('%Y%m%d_%H%M%S'))
            work.mkdir()
            expected=set()
            raw1,raw2=work/'input_R1.fastq',work/'input_R2.fastq'
            n=0
            with open(raw1,'xb') as x, open(raw2,'xb') as y:
                while pair is not None:
                    key=name(pair[0])
                    if key in expected:
                        raise ValueError(f'Duplicate input ID in {chunk}: {key!r}')
                    expected.add(key)
                    x.writelines(pair[0]); y.writelines(pair[1])
                    n+=1
                    if n>=a.chunk_pairs:
                        break
                    pair=paired(r1,r2)
            print(f'RUN: {chunk}, {n:,} pairs',flush=True)
            prefixes={label:work/label for label in ('rrna','non_rrna')}
            command=[binary,'--ref',str(ref),'--reads',str(raw1),'--reads',str(raw2),
                     '--workdir',str(work/'smr_work'),'--aligned',str(prefixes['rrna']),
                     '--other',str(prefixes['non_rrna']),'--fastx','--paired_in','--out2',
                     '--threads',str(a.threads)]
            (report/f'{chunk}.command.json').write_text(json.dumps(command,indent=2)+'\n')
            with open(report/f'{chunk}.console.log','w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            groups={label:(output_file(prefix,'fwd'),output_file(prefix,'rev'))
                    for label,prefix in prefixes.items()}
            counts=validate_groups(expected,groups)
            if sum(counts.values())!=n:
                raise ValueError(f'Count conservation failed for {chunk}')
            records=[]
            for label,files in groups.items():
                for mate,source in enumerate(files,1):
                    target=data/label/f'{chunk}_R{mate}.fastq.gz'
                    compress(source,target)
                    records.append(dict(group=label,mate=mate,path=str(target),
                                        size=target.stat().st_size,sha256=sha256(target)))
            for log in work.glob('*.log'):
                shutil.copy2(log,report/f'{chunk}.{log.name}')
            saved=dict(chunk=chunk,input_pairs=n,rrna_pairs=counts['rrna'],
                       non_rrna_pairs=counts['non_rrna'],files=records)
            atomic_json(done,saved)  # Only checkpoint after pairing checks and compressed file hashes.
            completed.append(saved); total_input+=n
            print(f'PASS: {chunk}: rRNA {counts["rrna"]:,}, non-rRNA {counts["non_rrna"]:,}',flush=True)
            # Only remove newly generated, validated scratch intermediates for THIS chunk.
            if work.parent!=scratch or not work.name.startswith(chunk+'.'):
                raise RuntimeError('Scratch cleanup path check failed')
            shutil.rmtree(work)
    if total_input!=a.expected_pairs:
        raise ValueError(f'Input count {total_input} differs from expected {a.expected_pairs}; outputs retained for review')
    total_rrna=sum(x['rrna_pairs'] for x in completed)
    total_non=sum(x['non_rrna_pairs'] for x in completed)
    if total_rrna+total_non!=total_input:
        raise ValueError('Global count mismatch')
    with open(report/'chunk_counts.tsv','w') as f:
        f.write('chunk\tinput_pairs\trrna_pairs\tnon_rrna_pairs\n')
        for c in completed:
            f.write(f'{c["chunk"]}\t{c["input_pairs"]}\t{c["rrna_pairs"]}\t{c["non_rrna_pairs"]}\n')
    for mate in (1,2):
        paths=[item['path'] for c in completed if c['non_rrna_pairs']
               for item in c['files'] if item['group']=='non_rrna' and item['mate']==mate]
        (report/f'non_rrna_R{mate}.list').write_text('\n'.join(paths)+'\n')
    atomic_json(report/'summary.json',dict(input_pairs=total_input,rrna_pairs=total_rrna,
                                          non_rrna_pairs=total_non,chunks=len(completed),
                                          annotation='rRNA-screened candidates; cross-species reference; residual rRNA/other origins require assessment'))
    with open(report/'SHA256SUMS','w') as f:
        for c in completed:
            for item in c['files']:
                f.write(f'{item["sha256"]}  {item["path"]}\n')
    print(f'COMPLETE: input={total_input:,}; rRNA={total_rrna:,} ({total_rrna/total_input:.2%}); '
          f'non-rRNA={total_non:,} ({total_non/total_input:.2%})',flush=True)
    print('Outputs:',data,'\nReports:',report,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/work/cyu/annotation/scotti-annotation'))
    p.add_argument('--scratch',type=Path,default=Path('/scratch/cyu/scotti_sortmerna_full_v1'))
    p.add_argument('--chunk-pairs',type=int,default=10_000_000)
    p.add_argument('--threads',type=int,default=12)
    p.add_argument('--expected-pairs',type=int,default=505_478_106)
    p.add_argument('--pilot-only',action='store_true')
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    try:
        main(args)
    except Exception as error:
        print(f'STOPPED: {error}\nOriginal inputs and completed outputs preserved. Review logs before resuming.',file=sys.stderr)
        raise SystemExit(1)
