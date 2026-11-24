#!/usr/bin/env python3
"""Freeze a verified working snapshot without changing any source annotation.

Standard library only. GenomeTools must be available. Original protein IDs and
coordinates are retained; the genome is recorded by path and SHA256, not copied.
"""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta(path, sequences=True):
    result = {}
    name = None
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                require(name not in result, 'Duplicate FASTA ID: '+name)
                result[name] = '' if sequences else 0
            elif line.strip():
                require(name is not None, 'Sequence before FASTA header')
                s = ''.join(line.split())
                if sequences:
                    result[name] += s
                else:
                    result[name] += len(s)
    require(bool(result) and all(result.values()), 'Empty FASTA sequence: '+str(path))
    return result


def gff(path, sizes):
    rows = []
    ids = {}
    with path.open() as handle:
        require(handle.readline().strip() == '##gff-version 3', 'Missing GFF3 header: '+str(path))
        for line in handle:
            if line.startswith('##FASTA'):
                raise ValueError('Embedded FASTA is not supported')
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\r\n').split('\t')
            require(len(f) == 9, 'Malformed GFF3 feature')
            require(f[0] in sizes and 1 <= int(f[3]) <= int(f[4]) <= sizes[f[0]],
                    'GFF3/genome coordinate mismatch')
            at = {k: v for item in f[8].split(';') if '=' in item
                  for k, v in [item.split('=', 1)]}
            if 'ID' in at:
                fid = unquote(at['ID'])
                if fid in ids:
                    old = ids[fid]
                    require(f[2] not in {'gene','mRNA'} and
                            (old[0][0:3], old[0][6], old[1].get('Parent')) ==
                            (f[0:3], f[6], at.get('Parent')),
                            'Inconsistent multipart GFF3 ID: '+fid)
                    require((int(f[3]),int(f[4])) not in old[2], 'Duplicate GFF3 segment: '+fid)
                    old[2].add((int(f[3]),int(f[4])))
                else:
                    ids[fid] = (f, at, {(int(f[3]),int(f[4]))})
            rows.append((f, at))
    for f, at in rows:
        if 'Parent' in at:
            require(all(unquote(v) in ids for v in at['Parent'].split(',')), 'Unresolved GFF3 Parent')
    return rows


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def check_proteins(protein_gff, all_proteins, coding, primary):
    genes = {unquote(at['ID']) for f, at in protein_gff if f[2] == 'gene'}
    tx = {unquote(at['ID']): unquote(at['Parent']) for f, at in protein_gff if f[2] == 'mRNA'}
    require(len(genes) == 9234 and len(tx) == 11464, 'Unexpected filtered7 gene/transcript count')
    require(set(tx.values()) == genes, 'Transcript/gene mismatch')
    require(set(all_proteins) == set(coding) == set(tx), 'GFF3/protein/CDS ID mismatch')
    require(len(primary) == 9234 and set(primary) <= set(tx), 'Unexpected primary protein set')
    require(len({tx[pid] for pid in primary}) == len(genes), 'Primary FASTA is not one protein per gene')
    require(all(primary[pid] == all_proteins[pid] for pid in primary), 'Primary protein sequence mismatch')
    require(not genes.intersection({'g1218','g196','g3336','g5529','g5573','g7051','g7437'}),
            'An excluded short gene is present')
    require(all('*' not in seq.rstrip('*') for seq in all_proteins.values()), 'Internal protein stop found')
    by_gene = defaultdict(list)
    for pid in tx:
        by_gene[tx[pid]].append(pid)
    for pid in primary:
        require(len(primary[pid].rstrip('*')) == max(len(all_proteins[p].rstrip('*')) for p in by_gene[tx[pid]]),
                'Selected primary protein is not a longest isoform: '+pid)
    return genes, tx


def busco_states(path):
    states = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\r\n').split('\t')
            require(len(f) >= 2 and f[1] in {'Complete','Duplicated','Fragmented','Missing'},
                    'Unexpected BUSCO full_table record')
            require(f[0] not in states or states[f[0]] == f[1], 'Conflicting BUSCO statuses')
            states[f[0]] = f[1]
    require(len(states) == 672, 'Expected metazoa_odb12 with 672 BUSCO groups')
    return states


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('/work/cyu/annotation/scotti-annotation'))
    p.add_argument('--genome', type=Path, default=Path('/work/cyu/assembly/worm/scotti/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa'))
    p.add_argument('--gt', type=Path, default=Path('/home/cyu/.conda/envs/genome_assembly/bin/gt'))
    p.add_argument('--out', type=Path)
    a = p.parse_args()
    result = a.root/'results'
    candidate = result/'filtered_candidates/Scotti_filtered7_v1'
    te = result/'annotation_qc/08_TE_protein_review_v1/review100_v2'
    ncrna = result/'ncrna_v1'
    out = a.out or result/'working_snapshots/Scotti_annotation_working_v1'
    partial = out.with_name(out.name+'.building')
    require(not out.exists() and not partial.exists(), 'Snapshot or staging directory exists; inspect it and use a new --out')
    require(a.gt.is_file(), 'GenomeTools executable missing: '+str(a.gt))

    protein_path = candidate/'Scotti_nuclear.filtered7.proteins.fa'
    cds_path = candidate/'Scotti_nuclear.filtered7.cds.fa'
    primary_path = candidate/'Scotti_nuclear.filtered7.longest.proteins.fa'
    gff_path = candidate/'Scotti_nuclear.filtered7.candidate.gff3'
    sensitivity_path = te/'Scotti_primary.without32_strong_TE_candidates.sensitivity.faa'
    drop_path = te/'Scotti_primary.strong32_TE_candidates.faa'
    review_path = te/'Scotti_TE.primary_case_review100.v2.tsv'
    nc_path = ncrna/'05_ncrna_candidate_v1/Scotti_nuclear.ncRNA.combined.candidate.gff3'
    context_path = ncrna/'06_ncrna_repeat_context_v1/summary.json'
    for path in [a.genome, protein_path, cds_path, primary_path, gff_path,
                 sensitivity_path, drop_path, review_path, nc_path, context_path]:
        require(path.is_file() and path.stat().st_size, 'Missing input: '+str(path))
    print('Checking genome, annotation IDs and protein versions...', flush=True)
    genome_hash = digest(a.genome)
    validated_paths = [protein_path, cds_path, primary_path, gff_path, sensitivity_path,
                       drop_path, review_path, nc_path, context_path]
    validated_hashes = {path: digest(path) for path in validated_paths}
    sizes = fasta(a.genome, sequences=False)
    require(len(sizes) == 412 and sum(sizes.values()) == 991317046 and 'Scotti_mt' not in sizes,
            'Unexpected selected RagTag nuclear genome')
    all_proteins, coding, primary = map(fasta, [protein_path, cds_path, primary_path])
    genes, tx = check_proteins(gff(gff_path, sizes), all_proteins, coding, primary)
    sensitivity, dropped = map(fasta, [sensitivity_path, drop_path])
    review = table(review_path)
    counts = Counter(r['review_group_v2'] for r in review)
    require(len(review) == 100 and len({r['protein_id'] for r in review}) == 100,
            'Expected 100 unique primary-protein TE cases')
    require(counts == {'strong_TE_candidate':32,'retain_TE_not_supported':21,'unresolved_keep_in_comparison':47},
            'Unexpected TE case counts')
    require(all(r['protein_id'] in primary and tx[r['protein_id']] == r['gene_id'] for r in review),
            'TE review is not linked to the current primary proteins')
    strong = {r['protein_id'] for r in review if r['review_group_v2'] == 'strong_TE_candidate'}
    require(set(dropped) == strong and len(sensitivity) == 9202 and
            set(sensitivity) == set(primary)-strong, 'Sensitivity subset does not omit exactly the reviewed 32 proteins')
    require(all(primary[pid] == seq for data in [sensitivity, dropped] for pid, seq in data.items()),
            'Sensitivity protein sequence mismatch')

    functional = candidate/'functional_evidence_v1'
    for name, expected in [('Scotti_filtered7.all_isoforms.evidence.tsv', all_proteins),
                           ('Scotti_filtered7.primary.evidence.tsv', primary)]:
        path = functional/name
        validated_hashes[path] = digest(path)
        evidence = table(path)
        require(len(evidence) == len(expected) and {r['protein_id'] for r in evidence} == set(expected),
                'Functional evidence/protein ID mismatch: '+name)
        require(all(tx[r['protein_id']] == r['gene_id'] for r in evidence),
                'Functional evidence/gene ID mismatch: '+name)

    nc = gff(nc_path, sizes)
    context = json.loads(context_path.read_text())
    nc_counts = Counter(f[2] for f, at in nc)
    require(len(nc) == 3681 and nc_counts['tRNA'] == 344, 'Unexpected ncRNA candidate counts')
    require(context['candidate_features'] == len(nc) and context['candidate_features_removed'] == 0,
            'ncRNA repeat audit count mismatch')

    baseline_table = candidate/'busco_longest_odb12/run_metazoa_odb12/full_table.tsv'
    sensitivity_table = te/'busco_without32_odb12/run_metazoa_odb12/full_table.tsv'
    for path in [baseline_table, sensitivity_table]:
        validated_hashes[path] = digest(path)
    before, after = map(busco_states, [baseline_table, sensitivity_table])
    require(set(before) == set(after), 'Different BUSCO group IDs between protein versions')
    changes = [dict(busco_id=bid, baseline=before[bid], without32=after[bid])
               for bid in sorted(before) if before[bid] != after[bid]]

    # Copy documented outputs only; no BAMs, raw reads, genome copies or search temp directories.
    plan = []
    def add_file(src, dest):
        require(src.is_file(), 'Missing snapshot source: '+str(src))
        plan.append((src, Path(dest)))
    for name in ['Scotti_nuclear.filtered7.candidate.gff3','Scotti_nuclear.filtered7.proteins.fa',
                 'Scotti_nuclear.filtered7.cds.fa','Scotti_nuclear.filtered7.longest.proteins.fa',
                 'excluded_models.tsv','excluded.proteins.fa','longest_isoform_selection.tsv','summary.json','README.txt']:
        add_file(candidate/name, Path('protein_annotation')/name)
    for name in ['Scotti_filtered7.all_isoforms.evidence.tsv','Scotti_filtered7.primary.evidence.tsv',
                 'Scotti_filtered7.swissprot.best_hit.tsv','Scotti_filtered7.swissprot.top5.tsv',
                 'Scotti_filtered7.emapper.annotations','Scotti_filtered7.interproscan.tsv',
                 'coverage_summary.json','README.txt','integrate_filtered7_evidence.py']:
        add_file(functional/name, Path('functional_evidence')/name)
    for name in ['Scotti_primary.without32_strong_TE_candidates.sensitivity.faa',
                 'Scotti_primary.strong32_TE_candidates.faa','Scotti_TE.primary_case_review100.v2.tsv',
                 'Scotti_primary.all_review_status.v2.tsv','strong32.all_isoform_context.tsv',
                 'summary.json','README.txt','finalize_scotti_te_case_review.py']:
        add_file(te/name, Path('TE_sensitivity')/name)
    for dirname in ['04_trna_candidate_v1','05_ncrna_candidate_v1','06_ncrna_repeat_context_v1']:
        folder = ncrna/dirname
        require(folder.is_dir(), 'Missing ncRNA source folder: '+str(folder))
        for src in sorted(folder.iterdir()):
            if src.is_file():
                add_file(src, Path('ncRNA')/dirname/src.name)
    for label, folder in [('baseline',candidate/'busco_longest_odb12'),('without32',te/'busco_without32_odb12')]:
        for src in sorted(folder.glob('short_summary*.txt')):
            add_file(src, Path('qc/BUSCO_odb12')/label/src.name)
        add_file(folder/'run_metazoa_odb12/full_table.tsv', Path('qc/BUSCO_odb12')/label/'full_table.tsv')
    require(len({dest for src, dest in plan}) == len(plan), 'Duplicate snapshot destination')
    partial.mkdir(parents=True)
    source_rows = []
    print(f'Copying and verifying {len(plan)} files...', flush=True)
    for src, dest in plan:
        dst = partial/dest
        dst.parent.mkdir(parents=True, exist_ok=True)
        initial = digest(src)
        require(src not in validated_hashes or initial == validated_hashes[src],
                'Validated source changed before copying: '+str(src))
        shutil.copy2(src, dst)
        require(digest(dst) == initial == digest(src), 'Source changed or copy failed: '+str(src))
        source_rows.append(dict(snapshot_path=str(dest), source_path=str(src), SHA256=initial,
                                bytes=dst.stat().st_size))
    qc = partial/'qc'
    qc.mkdir(exist_ok=True)
    for name, dest in [('protein',Path('protein_annotation')/gff_path.name),
                       ('ncRNA',Path('ncRNA/05_ncrna_candidate_v1')/nc_path.name)]:
        log = qc/(name+'.genometools_validation.log')
        with log.open('w') as handle:
            proc = subprocess.run([str(a.gt),'gff3validator','-typecheck','so',str(partial/dest)],
                                  stdout=handle, stderr=subprocess.STDOUT)
        require(proc.returncode == 0 and 'input is valid GFF3' in log.read_text(),
                'GenomeTools validation failed; inspect '+str(log))
    def tsv(path, fields, rows):
        with path.open('w') as handle:
            w = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            w.writeheader(); w.writerows(rows)
    tsv(qc/'BUSCO_status_changes.tsv',['busco_id','baseline','without32'],changes)
    tsv(partial/'source_manifest.tsv', ['snapshot_path','source_path','SHA256','bytes'], source_rows)
    require(genome_hash == digest(a.genome), 'Genome changed during snapshot construction')
    tsv(partial/'genome_reference.tsv',['source_path','SHA256','sequences','bp','included_in_snapshot'],
        [dict(source_path=str(a.genome), SHA256=genome_hash, sequences=len(sizes),
              bp=sum(sizes.values()), included_in_snapshot=False)])
    summary = dict(status='frozen_working_snapshot_not_final_annotation',
                   created_UTC=datetime.now(timezone.utc).isoformat(),
                   genes=len(genes), transcripts=len(tx), all_isoform_proteins=len(all_proteins),
                   primary_proteins=len(primary), sensitivity_primary_proteins=len(sensitivity),
                   strong_TE_candidates_omitted_in_sensitivity_only=len(strong),
                   unresolved_TE_case_models_retained=47,
                   ncRNA_candidates=len(nc), tRNA_candidates=344, Rfam_candidates=len(nc)-344,
                   ncRNA_TE_class_ge80pct=context['TE_class_ge80pct'],
                   BUSCO_baseline_status_counts=dict(Counter(before.values())),
                   BUSCO_without32_status_counts=dict(Counter(after.values())),
                   BUSCO_group_status_changes=len(changes),
                   protein_GFF3_validation='passed', ncRNA_GFF3_validation='passed',
                   source_files_modified=False, genome_copied=False,
                   comparative_input_replaced=False)
    (partial/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (partial/'STATUS.txt').write_text(
        'FROZEN WORKING SNAPSHOT; NOT A FINAL ANNOTATION RELEASE\n'
        'Baseline: filtered7, 9234 genes / 11464 isoforms / 9234 primary proteins.\n'
        'TE sensitivity: 9202 primary proteins, omitting only 32 reviewed strong TE-associated candidates.\n'
        'The 47 unresolved TE cases remain retained. Other primary proteins are not certified non-TE.\n'
        'ncRNA: 344 high-confidence tRNA loci and 3337 Rfam homology candidates; functional status is not universally validated.\n')
    (partial/'README.md').write_text(
        '# S. cotti annotation working snapshot\n\n'
        'This package freezes the current candidate annotation and its evidence. Original gene/protein IDs and coordinates are preserved.\n\n'
        '- `protein_annotation/`: filtered7 structural baseline, all protein/CDS isoforms, longest representatives and seven-model exclusion record.\n'
        '- `functional_evidence/`: Swiss-Prot, eggNOG and InterPro evidence linked to the baseline IDs. Evidence is not a validated product assignment.\n'
        '- `TE_sensitivity/`: reviewed case decisions and the 9202-record sensitivity FASTA. It is not a complete non-TE proteome or a replacement structural annotation.\n'
        '- `ncRNA/`: separate candidate GFF3 components, provenance and repeat-context audit. Protein and ncRNA GFF3 files are not merged here.\n'
        '- `qc/`: fresh GenomeTools checks and existing BUSCO full tables/summaries; BUSCO was not rerun.\n'
        '- `genome_reference.tsv`: exact nuclear genome path and SHA256. The genome itself is not copied into this package.\n'
        '- `source_manifest.tsv`: copied-file provenance; `SHA256SUMS` covers the snapshot files.\n\n'
        'The original source scripts and their README files retain their historical candidate-stage wording. Current snapshot checks are recorded in the top-level summary and qc directory. Nested SHA256SUMS/input.SHA256SUMS are historical provenance; use the top-level SHA256SUMS to verify this snapshot.\n\n'
        'Scope: no new gene prediction, protein search, TE-domain search, RNA alignment, BUSCO run or comparative OrthoFinder run was performed. Comparative inputs and existing results remain unchanged.\n\n'
        'Known review items: 47 unresolved TE cases; additional short models; possible missed/repaired genes; functional ncRNA versus pseudogene status, including high-copy U1. Regional RNA coverage alone does not validate a protein model.\n\n'
        'For baseline comparative analysis use `protein_annotation/Scotti_nuclear.filtered7.longest.proteins.fa`; for the specific 32-candidate sensitivity comparison use the clearly named FASTA under `TE_sensitivity/`. These are distinct inputs.\n\n'
        'Verify from this directory with `sha256sum -c SHA256SUMS`.\n')
    shutil.copy2(Path(__file__), partial/Path(__file__).name)
    with (partial/'SHA256SUMS').open('w') as handle:
        for path in sorted(partial.rglob('*')):
            if path.is_file() and path != partial/'SHA256SUMS':
                handle.write(f'{digest(path)}  {path.relative_to(partial)}\n')
    # Publish the completed directory only after all checks pass.
    require(not out.exists(), 'Snapshot destination appeared during construction; staging preserved')
    partial.rename(out)
    print('Snapshot:', out)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
