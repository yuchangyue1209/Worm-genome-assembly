#!/usr/bin/env python3
"""Exercise publication/copy/checksum/permission logic with synthetic data.

The external GenomeTools validator and nuclear size check are fixture stubs.
Real validation and real genome size checks remain mandatory in the server script.
"""
import contextlib
import grp
import io
import json
import os
import stat
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import publish_scotti_shared_release as m


def put(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def manifest(root):
    with (root/'SHA256SUMS').open('w') as handle:
        for path in sorted(root.rglob('*')):
            if path.is_file() and path != root/'SHA256SUMS':
                handle.write(m.checksum_line(path, root))


def verify(root, log):
    for raw in (root/'SHA256SUMS').read_text().splitlines():
        escaped = raw.startswith('\\')
        raw = raw[1:] if escaped else raw
        expected, name = raw.split('  ',1)
        if escaped:
            decoded = []
            i = 0
            while i < len(name):
                if name[i] == '\\':
                    i += 1
                    decoded.append('\n' if name[i] == 'n' else name[i])
                else:
                    decoded.append(name[i])
                i += 1
            name = ''.join(decoded)
        assert m.sha(root/name) == expected, name
    log.write_text('Fixture checksum verification passed\n')


def validator(gt, path, log):
    assert path.read_text().startswith('##gff-version 3\n')
    log.write_text('input is valid GFF3\n')


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    src = root/'snapshot'
    mt = root/'mt'
    shared = root/'shared'
    shared.mkdir()
    genome = root/'genome.fa'
    put(genome, '>chr\nACGT\n')
    summary = dict(status='frozen_working_snapshot_not_final_annotation', genes=9234,
                   transcripts=11464, all_isoform_proteins=11464, primary_proteins=9234, sensitivity_primary_proteins=9202,
                   ncRNA_candidates=3681, BUSCO_group_status_changes=0,
                   protein_GFF3_validation='passed', ncRNA_GFF3_validation='passed')
    put(src/'summary.json', json.dumps(summary))
    put(src/'genome_reference.tsv', 'source_path\tSHA256\n'+str(genome)+'\t'+m.sha(genome)+'\n')
    for name in ['Scotti_nuclear.filtered7.proteins.fa','Scotti_nuclear.filtered7.longest.proteins.fa',
                 'Scotti_nuclear.filtered7.cds.fa']:
        put(src/'protein_annotation'/name, '>g1.t1\nMAA\n')
    put(src/'protein_annotation/Scotti_nuclear.filtered7.candidate.gff3', '##gff-version 3\n')
    put(src/'ncRNA/05_ncrna_candidate_v1/Scotti_nuclear.ncRNA.combined.candidate.gff3', '##gff-version 3\n')
    put(src/'TE_sensitivity/Scotti_primary.without32_strong_TE_candidates.sensitivity.faa', '>sensitivity\nMAA\n')
    put(src/'functional_evidence/Scotti_filtered7.primary.evidence.tsv','protein_id\tgene_id\ng1.t1\tg1\n')
    put(src/'README.txt', 'Frozen working source\n')
    manifest(src)
    put(mt/'Scotti_mitochondrial.fa', '>contig_1\n'+'A'*15437+'\n')
    gff = '##gff-version 3\n'
    for i in range(36):
        root_kind = 'gene' if i < 12 else 'ncRNA_gene'
        child_kind = 'mRNA' if i < 12 else ('tRNA' if i < 34 else 'rRNA')
        gff += f'contig_1\tcurated\t{root_kind}\t{i*10+1}\t{i*10+9}\t.\t+\t.\tID=gene{i}\n'
        gff += f'contig_1\tcurated\t{child_kind}\t{i*10+1}\t{i*10+9}\t.\t+\t.\tID=tx{i};Parent=gene{i}\n'
        if i < 12:
            gff += f'contig_1\tcurated\tCDS\t{i*10+1}\t{i*10+9}\t.\t+\t0\tID=CDS{i};Parent=tx{i}\n'
    put(mt/'Scotti_mitochondrial.curated_candidate.gff3', gff)
    counts = m.mitochondrial_counts(mt/'Scotti_mitochondrial.curated_candidate.gff3')
    assert counts['gene']==12 and counts['ncRNA_gene']==24
    put(mt/'review\\notes.tsv','test\n')  # Check GNU manifest escaping.
    manifest(mt)
    gt = root/'gt_stub'
    put(gt,'fixture placeholder\n')
    original = {p:m.sha(p) for folder in [src,mt] for p in folder.rglob('*') if p.is_file()}
    dest = shared/'release'
    real_size = m.fasta_size
    def sizes(path):
        return (412,991317046) if path.name=='Scotti_nuclear.fa' else real_size(path)
    args = [str(m.__file__),'--snapshot',str(src),'--mitochondria',str(mt),
            '--destination',str(dest),'--gt',str(gt),'--group',grp.getgrgid(os.getgid()).gr_name]
    with patch.object(sys,'argv',args), patch.object(m,'verify_manifest',verify), \
         patch.object(m,'validate_gff',validator), patch.object(m,'fasta_size',sizes), \
         patch.object(m.shutil,'which',return_value='fixture_stub'), contextlib.redirect_stdout(io.StringIO()):
        m.main()
        assert dest.is_dir() and not dest.with_name(dest.name+'.building').exists()
        result = json.loads((dest/'summary.json').read_text())
        assert result['status']=='scientific_release_v1.0_preNCBI'
        assert result['primary_proteins']==9234 and result['all_isoform_proteins']==11464
        assert len(list((dest/'proteins').glob('*.fa')))==2
        assert (dest/'working_snapshot/TE_sensitivity/Scotti_primary.without32_strong_TE_candidates.sensitivity.faa').is_file()
        verify(dest,shared/'fixture_verify.log')
        assert all(m.sha(p)==h for p,h in original.items())
        published = dest/'proteins/Scotti_proteins.longest_isoform.v1.0.fa'
        source = src/'protein_annotation/Scotti_nuclear.filtered7.longest.proteins.fa'
        assert published.stat().st_ino != source.stat().st_ino
        assert stat.S_IMODE(published.stat().st_mode) & stat.S_IRGRP
        if sys.platform != 'darwin':
            assert stat.S_IMODE(dest.stat().st_mode) & stat.S_ISGID
        try:
            m.main()
        except ValueError as error:
            assert 'Destination/staging exists' in str(error)
        else:
            raise AssertionError('Existing release was not protected')
    # Simulate an interrupted publication with the existing copied genome/packages.
    staging = dest.with_name(dest.name+'.building')
    dest.rename(staging)
    nuclear_inode = (staging/'assembly/Scotti_nuclear.fa').stat().st_ino
    with patch.object(sys,'argv',args+['--resume']), patch.object(m,'verify_manifest',verify), \
         patch.object(m,'validate_gff',validator), patch.object(m,'fasta_size',sizes), \
         patch.object(m.shutil,'which',return_value='fixture_stub'), contextlib.redirect_stdout(io.StringIO()):
        m.main()
    assert (dest/'assembly/Scotti_nuclear.fa').stat().st_ino==nuclear_inode
    assert all(m.sha(p)==h for p,h in original.items())
    verify(dest,shared/'fixture_resume_verify.log')
    mismatched = root/'mismatched_copy.fa'
    put(mismatched,'preserve this\n')
    try:
        m.copy_verified(genome,mismatched)
    except ValueError as error:
        assert 'differs from source' in str(error)
    else:
        raise AssertionError('A mismatching existing staged file was overwritten')
    assert mismatched.read_text()=='preserve this\n'
print('PASS: publication layout, selected baseline, sensitivity kept only in provenance,')
print('      copy/checksum integrity, escaped filenames, shared permissions and overwrite protection.')
print('PASS: 12 gene + 24 ncRNA_gene mitochondrial schema, verified resume without recopying,')
print('      and refusal to overwrite mismatching staged data.')
print('NOTE: real GenomeTools validation and production genome counts run on the server.')
