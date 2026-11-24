#!/usr/bin/env python3
"""Publish the selected S. cotti baseline as v1.0 pre-NCBI in /work/shared.

Original files are preserved. Existing destinations are never overwritten.
The exact nuclear genome and the mitochondrial working package are included.
"""
import argparse
import csv
import grp
import hashlib
import json
import os
import shutil
import stat
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def checksum_line(path, base):
    name = str(path.relative_to(base))
    escaped = '\\' in name or '\n' in name
    if escaped:
        name = name.replace('\\', '\\\\').replace('\n', '\\n')
    return ('\\' if escaped else '')+sha(path)+'  '+name+'\n'


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def tsv(path, fields, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def copy_verified(src, dst):
    require(src.is_file() and not src.is_symlink(), 'Missing source or unsupported symlink: '+str(src))
    expected = sha(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if os.path.lexists(dst):
        require(dst.is_file() and not dst.is_symlink() and sha(dst) == expected == sha(src),
                'Existing staged file differs from source; preserved for inspection: '+str(dst))
        return expected
    shutil.copy2(src, dst)
    require(sha(dst) == expected == sha(src), 'Copy mismatch/source changed: '+str(src))
    return expected


def copy_tree(src, dst):
    require(src.is_dir() and not src.is_symlink(), 'Missing source directory: '+str(src))
    for path in sorted(src.rglob('*')):
        require(not path.is_symlink(), 'Source symlink requires explicit handling: '+str(path))
        if path.is_file():
            copy_verified(path, dst/path.relative_to(src))
        elif path.is_dir():
            (dst/path.relative_to(src)).mkdir(parents=True, exist_ok=True)


def verify_manifest(root, log):
    require((root/'SHA256SUMS').is_file(), 'Checksum manifest missing: '+str(root))
    with log.open('w') as handle:
        proc = subprocess.run(['sha256sum','-c','SHA256SUMS'], cwd=root,
                              stdout=handle, stderr=subprocess.STDOUT)
    require(proc.returncode == 0, 'Checksum verification failed; inspect '+str(log))


def alias(src, dst):
    """Alias only newly copied files inside this shared package, never originals."""
    require(src.is_file(), 'Alias source missing: '+str(src))
    if os.path.lexists(dst):
        require(dst.is_file() and not dst.is_symlink() and sha(dst) == sha(src),
                'Existing alias differs from copied source; preserved: '+str(dst))
        return 'internal_hardlink' if os.path.samefile(src, dst) else 'verified_copy'
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)
        method = 'internal_hardlink'
    except OSError:
        copy_verified(src, dst)
        method = 'verified_copy'
    return method


def fasta_size(path):
    records = 0
    bases = 0
    name = None
    ids = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                require(name not in ids, 'Duplicate FASTA ID')
                ids.add(name); records += 1
            elif line.strip():
                require(name is not None, 'Sequence before FASTA header')
                bases += len(''.join(line.split()))
    return records, bases


def validate_gff(gt, path, log):
    with log.open('w') as handle:
        proc = subprocess.run([str(gt),'gff3validator','-typecheck','so',str(path)],
                              stdout=handle, stderr=subprocess.STDOUT)
    require(proc.returncode == 0 and 'input is valid GFF3' in log.read_text(),
            'GFF3 validation failed; inspect '+str(log))


def mitochondrial_counts(path):
    """Count top-level loci across gene AND ncRNA_gene, not transcript rows."""
    counts = Counter()
    ids = set()
    for line in path.read_text().splitlines():
        if line.startswith('##FASTA'):
            break
        if line.startswith('#') or not line.strip():
            continue
        f = line.split('\t')
        require(len(f) == 9, 'Malformed mitochondrial GFF3')
        counts[f[2]] += 1
        if f[2] in {'gene','ncRNA_gene'}:
            attrs = dict(item.split('=',1) for item in f[8].split(';') if '=' in item)
            fid = attrs.get('ID')
            require(fid and fid not in ids and 'Parent' not in attrs,
                    'Missing/duplicate or non-root mitochondrial gene ID')
            ids.add(fid)
    require(len(ids) == counts['gene']+counts['ncRNA_gene'] == 36 and
            counts['mRNA'] == 12 and counts['CDS'] == 12 and
            counts['tRNA'] == 22 and counts['rRNA'] == 2,
            'Unexpected mitochondrial locus/feature counts: '+json.dumps(dict(counts),sort_keys=True))
    return dict(counts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot', type=Path, default=Path('/work/cyu/annotation/scotti-annotation/results/working_snapshots/Scotti_annotation_working_v1'))
    p.add_argument('--mitochondria', type=Path, default=Path('/work/cyu/assembly/worm/scotti/mitochondria/08_working_release/Scotti_mt_working_v1'))
    p.add_argument('--destination', type=Path, default=Path('/work/shared/Scotti_release_v1.0_preNCBI_20261118'))
    p.add_argument('--group', default='shared')
    p.add_argument('--gt', type=Path, default=Path('/home/cyu/.conda/envs/genome_assembly/bin/gt'))
    p.add_argument('--resume', action='store_true',
                   help='Verify/reuse the existing .building copy; mismatching staged data are not overwritten')
    a = p.parse_args()
    src = a.snapshot.resolve()
    mt = a.mitochondria.resolve()
    dest = a.destination.absolute()
    staging = dest.with_name(dest.name+'.building')
    verification_log = dest.with_name(dest.name+'.verification.log')
    require(not os.path.lexists(dest), 'Destination/staging exists; completed release will not be overwritten')
    if a.resume:
        require(staging.is_dir() and not staging.is_symlink(), '--resume requires the existing .building directory')
        require(not any(path.is_symlink() for path in staging.rglob('*')),
                'Existing staged symlink requires inspection before resuming')
    else:
        require(not os.path.lexists(staging), 'Staging exists; use --resume to verify and reuse it')
    require(not os.path.lexists(verification_log) or
            (a.resume and verification_log.is_file() and not verification_log.is_symlink()),
            'Verification log already exists; choose a new --destination')
    resolved_dest = dest.resolve()
    require(all(source != resolved_dest and source not in resolved_dest.parents and resolved_dest not in source.parents
                for source in [src,mt]), 'Source and destination must be separate directories')
    require(dest.parent.is_dir(), 'Shared parent directory missing: '+str(dest.parent))
    require(a.gt.is_file() and shutil.which('sha256sum'), 'GenomeTools or sha256sum missing')
    gid = grp.getgrnam(a.group).gr_gid
    summary = json.loads((src/'summary.json').read_text())
    require(summary['status'] == 'frozen_working_snapshot_not_final_annotation', 'Unexpected source status')
    expected = {'genes':9234,'transcripts':11464,'all_isoform_proteins':11464,'primary_proteins':9234,
                'sensitivity_primary_proteins':9202,'ncRNA_candidates':3681,
                'BUSCO_group_status_changes':0}
    require(all(summary.get(k) == value for k,value in expected.items()), 'Source snapshot count mismatch')
    require(summary['protein_GFF3_validation'] == summary['ncRNA_GFF3_validation'] == 'passed',
            'Source GFF3 validation status mismatch')
    genome_rows = table(src/'genome_reference.tsv')
    require(len(genome_rows) == 1, 'Expected one nuclear genome reference')
    genome = Path(genome_rows[0]['source_path'])
    require(genome.is_file() and sha(genome) == genome_rows[0]['SHA256'], 'Nuclear genome/reference checksum mismatch')
    mt_fa = mt/'Scotti_mitochondrial.fa'
    mt_gff = mt/'Scotti_mitochondrial.curated_candidate.gff3'
    require(mt_fa.is_file() and mt_gff.is_file() and (mt/'SHA256SUMS').is_file(),
            'Mitochondrial working package missing')
    require(fasta_size(mt_fa) == (1,15437), 'Unexpected mitochondrial FASTA length/count')
    mt_summary = mitochondrial_counts(mt_gff)

    staging.mkdir(exist_ok=a.resume)
    qc = staging/'quality_control'
    qc.mkdir(exist_ok=a.resume)
    print('Verifying the selected baseline snapshot and mitochondrial package...', flush=True)
    verify_manifest(src, qc/'source_snapshot.checksums.log')
    verify_manifest(mt, qc/'source_mitochondria.checksums.log')
    print(('Resuming: verifying and reusing copied files...' if a.resume else
           'Copying the snapshot and genome to shared staging...'), flush=True)
    embedded = staging/'working_snapshot'
    copy_tree(src, embedded)
    verify_manifest(embedded, qc/'copied_snapshot.checksums.log')
    nuclear_dest = staging/'assembly/Scotti_nuclear.fa'
    nuclear_hash = copy_verified(genome, nuclear_dest)
    require(nuclear_hash == genome_rows[0]['SHA256'], 'Copied nuclear genome/reference mismatch')
    require(fasta_size(nuclear_dest) == (412,991317046), 'Copied nuclear FASTA size mismatch')
    mt_dest = staging/'assembly/mitochondria_working_v1'
    copy_tree(mt, mt_dest)
    verify_manifest(mt_dest, qc/'copied_mitochondria.checksums.log')

    aliases = [
        (embedded/'protein_annotation/Scotti_nuclear.filtered7.candidate.gff3',
         Path('annotation/Scotti_nuclear.protein_coding.v1.0.gff3')),
        (embedded/'protein_annotation/Scotti_nuclear.filtered7.cds.fa',
         Path('annotation/Scotti_nuclear.CDS.all_isoforms.v1.0.fa')),
        (embedded/'ncRNA/05_ncrna_candidate_v1/Scotti_nuclear.ncRNA.combined.candidate.gff3',
         Path('annotation/Scotti_nuclear.ncRNA.candidate.v1.0.gff3')),
        (embedded/'protein_annotation/Scotti_nuclear.filtered7.proteins.fa',
         Path('proteins/Scotti_proteins.all_isoforms.v1.0.fa')),
        (embedded/'protein_annotation/Scotti_nuclear.filtered7.longest.proteins.fa',
         Path('proteins/Scotti_proteins.longest_isoform.v1.0.fa')),
        (mt_dest/mt_fa.name, Path('assembly/Scotti_mitochondrial.fa')),
        (mt_dest/mt_gff.name, Path('assembly/Scotti_mitochondrial.candidate.v1.0.gff3')),
    ]
    for path in sorted((embedded/'functional_evidence').iterdir()):
        if path.is_file():
            aliases.append((path,Path('functional_annotation')/path.name))
    alias_rows = []
    for source, target in aliases:
        method = alias(source, staging/target)
        alias_rows.append(dict(shared_path=str(target), copied_snapshot_path=str(source.relative_to(staging)),
                               method=method, SHA256=sha(staging/target)))
    tsv(staging/'file_map.tsv', ['shared_path','copied_snapshot_path','method','SHA256'], alias_rows)
    for label, path in [
        ('nuclear_protein',staging/'annotation/Scotti_nuclear.protein_coding.v1.0.gff3'),
        ('nuclear_ncRNA',staging/'annotation/Scotti_nuclear.ncRNA.candidate.v1.0.gff3'),
        ('mitochondrial',staging/'assembly/Scotti_mitochondrial.candidate.v1.0.gff3'),
    ]:
        validate_gff(a.gt,path,qc/(label+'.genometools_validation.log'))
    publication = dict(status='scientific_release_v1.0_preNCBI',
                       copied_UTC=datetime.now(timezone.utc).isoformat(),
                       destination=str(dest), source_snapshot=str(src), source_mitochondria=str(mt),
                       source_nuclear_genome=str(genome), nuclear_genome_SHA256=nuclear_hash,
                       genes=9234, all_isoform_proteins=11464, primary_proteins=9234,
                       selected_protein_baseline='filtered7_all_9234_genes_retained',
                       strong_TE_candidate_genes_retained_in_release=32,
                       sensitivity_primary_proteins_in_provenance_only=9202,
                       nuclear_ncRNA_candidates=3681,
                       mitochondrial_bp=15437, mitochondrial_genes=36,
                       mitochondrial_feature_counts=mt_summary,
                       shared_group=a.group, source_files_modified=False,
                       comparative_input_replaced=False)
    (staging/'summary.json').write_text(json.dumps(publication,indent=2)+'\n')
    (staging/'STATUS.txt').write_text(
        'SCIENTIFIC RELEASE v1.0; PRE-NCBI; NO ACCESSION ASSIGNED\n'
        'Selected release: 9234 genes, 11464 protein isoforms, 9234 longest representatives.\n'
        'All 32 strong TE-associated candidate genes remain in the selected baseline.\n'
        'The 9202-protein sensitivity set is retained only in working_snapshot/TE_sensitivity/ as provenance.\n'
        'Nuclear ncRNA features and the mitochondrial annotation retain their documented candidate status.\n')
    (staging/'README.md').write_text(
        '# Schistocephalus cotti release v1.0 (pre-NCBI)\n\n'
        'The user selected the filtered7 baseline for release: 9234 genes, 11464 protein isoforms and 9234 longest representatives. This follows the assembly/annotation/proteins/functional_annotation/quality_control layout used for the S. solidus shared release. This is a scientific data release; no NCBI accession is assigned. Original files and comparative-analysis inputs remain unchanged.\n\n'
        '- `assembly/`: the exact 412-sequence, 991317046-bp RagTag nuclear genome, plus the separate 15437-bp mitochondrial genome/annotation. The full mitochondrial working package is preserved in `assembly/mitochondria_working_v1/`.\n'
        '- `annotation/`: separate protein-coding GFF3, ncRNA candidate GFF3 and all-isoform CDS FASTA. No combined nuclear GFF3 or GTF has been synthesized.\n'
        '- `proteins/`: the selected release protein sets: all 11464 isoforms and 9234 longest representatives.\n'
        '- `functional_annotation/`: preserved Swiss-Prot/eggNOG/InterPro evidence. These are evidence tables, not a curated product-name release.\n'
        '- `quality_control/`: source/copied checksum checks and fresh validation of the three shared GFF3 files.\n'
        '- `working_snapshot/`: the exact frozen source snapshot with its original directory structure, TE case decisions, the 9202-record sensitivity set, ncRNA repeat-context audit, BUSCO tables and provenance. Its historical WORKING/CANDIDATE status files describe the source stage; the top-level STATUS.txt and summary.json identify this selected release.\n\n'
        'The standard-layout files are internal hardlinks to the newly copied packages when supported, with verified-copy fallback. They never link to the original /work/cyu files. The top-level file_map.tsv lists the mappings. Source paths embedded in historical manifests/logs are provenance; included data are available locally in this shared directory.\n\n'
        'The release retains all 9234 primary proteins, including the 32 strong TE-associated candidate genes and the 47 unresolved TE cases. The 9202-record sensitivity file is provenance only, not the selected release proteome. It omits only the reviewed 32 candidates and is not a complete non-TE proteome. The nuclear ncRNA catalogue is a separate candidate component: 344 high-confidence tRNA loci and 3337 Rfam homology candidates; their functional/pseudogene status is not universally resolved. The mitochondrial annotation is preserved with its original curated-candidate status.\n\n'
        'Source: /work/shared/Ssolidus_release_v1.0_preNCBI_20260831 was used as a layout reference only; no S. solidus data are copied into this package.\n\n'
        'Group permissions follow the solidus convention: shared group, group read/write, directory traversal and setgid on package directories.\n\n'
        'Verify this entire package with `sha256sum -c SHA256SUMS` from its root. The final verification log is written alongside this directory as <package>.verification.log.\n')
    shutil.copy2(Path(__file__), staging/Path(__file__).name)
    expected_files = {Path('working_snapshot')/path.relative_to(src)
                      for path in src.rglob('*') if path.is_file()}
    expected_files |= {Path('assembly/mitochondria_working_v1')/path.relative_to(mt)
                       for path in mt.rglob('*') if path.is_file()}
    expected_files |= {target for source,target in aliases}
    expected_files |= {Path('assembly/Scotti_nuclear.fa'),Path('file_map.tsv'),Path('summary.json'),
                       Path('STATUS.txt'),Path('README.md'),Path(__file__).name}
    expected_files = {Path(path) for path in expected_files}
    for name in ['source_snapshot.checksums.log','source_mitochondria.checksums.log',
                 'copied_snapshot.checksums.log','copied_mitochondria.checksums.log',
                 'nuclear_protein.genometools_validation.log','nuclear_ncRNA.genometools_validation.log',
                 'mitochondrial.genometools_validation.log']:
        expected_files.add(Path('quality_control')/name)
    unexpected = [str(path.relative_to(staging)) for path in staging.rglob('*')
                  if path.is_file() and path.relative_to(staging) not in expected_files | {Path('SHA256SUMS')}]
    require(not unexpected, 'Unexpected staged files preserved; inspect: '+', '.join(unexpected[:10]))
    require(not any(path.is_symlink() for path in staging.rglob('*')), 'Staged symlink requires review')
    with (staging/'SHA256SUMS').open('w') as handle:
        for path in sorted(staging.rglob('*')):
            if path.is_file() and path != staging/'SHA256SUMS':
                handle.write(checksum_line(path,staging))

    # Permission changes apply to the new shared copy only.
    for path in [staging]+sorted(staging.rglob('*')):
        os.chown(path,-1,gid)
        mode = stat.S_IMODE(path.stat().st_mode) | stat.S_IRGRP | stat.S_IWGRP
        if path.is_dir():
            mode |= stat.S_IXGRP | stat.S_ISGID
        elif mode & stat.S_IXUSR:
            mode |= stat.S_IXGRP
        os.chmod(path,mode)
    verify_manifest(staging,verification_log)
    os.chown(verification_log,-1,gid)
    os.chmod(verification_log,stat.S_IMODE(verification_log.stat().st_mode)|stat.S_IRGRP|stat.S_IWGRP)
    require(not os.path.lexists(dest), 'Destination appeared; staging preserved')
    staging.rename(dest)
    print('PASS: selected v1.0 pre-NCBI release created:',dest)
    print('PASS: source/copied checksums and all three GFF3 files validated.')
    print('PASS: shared group permissions applied; original sources preserved.')
    print(json.dumps(publication,indent=2))


if __name__ == '__main__':
    main()
