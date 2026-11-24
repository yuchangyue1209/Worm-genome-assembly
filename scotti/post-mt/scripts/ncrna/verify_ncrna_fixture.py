#!/usr/bin/env python3
"""Synthetic checks of candidate selection, overlap policy and coordinate handling."""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

script = Path(__file__).with_name('build_scotti_ncrna_candidate.py')
with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    base = root/'results/ncrna_v1'
    (base/'02_rfam').mkdir(parents=True)
    (base/'04_trna_candidate_v1').mkdir()
    pg = root/'results/filtered_candidates/Scotti_filtered7_v1'
    pg.mkdir(parents=True)
    genome = root/'genome.fa'
    genome.write_text('>chr\n'+'A'*3000+'\n')
    trna = base/'04_trna_candidate_v1/Scotti_nuclear.tRNA.high_confidence.candidate.gff3'
    trna.write_text('##gff-version 3\n'+''.join(
        f'chr\ttRNAscan-SE\ttRNA\t{s}\t{e}\t80\t{strand}\t.\tID=t{i}\n'
        for i, (s, e, strand) in enumerate([(10,30,'+'), (40,60,'-'), (70,90,'+')],1)))
    protein = pg/'Scotti_nuclear.filtered7.candidate.gff3'
    protein.write_text('##gff-version 3\n'
        'chr\tBRAKER\tCDS\t1000\t1080\t.\t+\t0\tParent=p1\n'
        'chr\tBRAKER\tCDS\t1050\t1090\t.\t+\t0\tParent=p2\n')
    lines = []
    def add(family, acc, s, e, strand='+', score=100, ev='1e-20', olp='*', mf=1, mt=40, ml=40, trunc='no'):
        a, b = (s,e) if strand=='+' else (e,s)
        fields = [str(len(lines)+1), family, acc, 'chr', '-', '-', 'cm',str(mf),str(mt),
                  str(a),str(b),strand,trunc,'1','0.5','0',str(score),ev,'!',olp,
                  '-','-','-','-','-','-',str(ml),'3000','test description']
        lines.append(' '.join(fields)+'\n')
    add('U2','RF00004',100,150)
    add('U2','RF00004',160,220,'-')
    add('U5','RF00020',105,155,score=80,ev='1e-10')
    add('SSU_rRNA_eukarya','RF01960',300,350,olp='=')
    add('U4','RF00015',400,420,mt=21)
    add('Plant_SRP','RF01855',500,550)
    add('mir-71','RF00832',600,650)
    add('tRNA-Sec','RF01852',700,740)
    add('U6','RF00026',20,40)
    add('U11','RF00548',800,850,olp='$')
    add('U1','RF00003',1000,1100,score=300,ev='0')
    add('U2','RF00004',1005,1105,score=250,ev='0')
    add('U4','RF00015',1200,1300,score=10,ev='1e-2000')
    add('U2','RF00004',1200,1300,score=100,ev='1e-1000')
    add('U2','RF00004',100,150,'-')
    (base/'02_rfam/Scotti_nuclear.rfam.tblout').write_text('# fixture\n'+''.join(lines))
    cmd = [sys.executable,str(script),'--root',str(root),'--genome',str(genome),'--expected-trna','3']
    subprocess.run(cmd,check=True,stdout=subprocess.PIPE,text=True)
    out = base/'05_ncrna_candidate_v1'
    summary = json.loads((out/'summary.json').read_text())
    assert summary['raw_Rfam_hits']==15
    assert summary['Rfam_core_candidate_features']==7
    assert summary['combined_candidate_features']==10
    with (out/'Rfam_all_hits.review.tsv').open() as f:
        audit = list(csv.DictReader(f,delimiter='\t'))
    assert audit[2]['reasons'].startswith('residual_same_strand_candidate_overlap:')
    assert 'lower_scoring_overlap' in audit[3]['reasons']
    assert 'short_model_span' in audit[4]['reasons']
    assert 'outside_project_core' in audit[5]['reasons']
    assert audit[7]['decision']=='crosscheck_only_tRNAscan_is_primary'
    assert 'overlaps_curated_tRNA' in audit[8]['reasons']
    assert audit[9]['decision']=='core_candidate'  # $ is not silently discarded.
    assert audit[10]['same_strand_CDS_overlap_bp']=='91'  # Union, not sum.
    assert audit[10]['decision']=='core_candidate' and audit[11]['decision']=='review_only'
    assert audit[12]['decision']=='core_candidate' and audit[13]['decision']=='review_only'
    assert audit[14]['decision']=='core_candidate'  # Opposite strand is independent.
    gff = (out/'Scotti_nuclear.ncRNA.combined.candidate.gff3').read_text()
    assert 'miRNA_precursor_homology_only' in gff
    bed = [line.split('\t') for line in (out/'Scotti_nuclear.ncRNA.combined.candidate.bed').read_text().splitlines()]
    assert any(f[1:3]==['159','220'] and f[5]=='-' for f in bed)
    # Audit classified TE context separately from direct small-RNA matches.
    rm = root/'repeats.out'
    rm.write_text('   SW perc div del ins query\n'+''.join(
        f'100 0.0 0.0 0.0 chr {s} {e} (0) + repeat{i} {kind} 1 50 (0) {i}\n'
        for i, (s,e,kind) in enumerate([
            (1000,1060,'LINE/CR1'), (1050,1080,'LTR/Gypsy'),
            (100,150,'snRNA'), (10,30,'tRNA'), (600,620,'Unknown'),
            (1200,1220,'Simple_repeat'), (1230,1240,'Low_complexity')],1)))
    audit_script = script.with_name('audit_scotti_ncrna_repeat_context.py')
    gff_before_audit = (out/'Scotti_nuclear.ncRNA.combined.candidate.gff3').read_bytes()
    subprocess.run([sys.executable,str(audit_script),'--root',str(root),
                    '--genome',str(genome),'--repeatmasker',str(rm)],
                   check=True,stdout=subprocess.PIPE,text=True)
    ao = base/'06_ncrna_repeat_context_v1'
    context = json.loads((ao/'summary.json').read_text())
    assert context['candidate_features']==10 and context['TE_class_ge80pct']==1
    assert context['small_RNA_ge80pct']==3 and context['candidate_features_removed']==0
    with (ao/'ncRNA_repeat_context.by_feature.tsv').open() as f:
        feature_context = list(csv.DictReader(f,delimiter='\t'))
    u1 = next(r for r in feature_context if r['family']=='U1')
    assert u1['TE_class_union_bp']=='81'  # Union, not 61+31.
    assert all(r['TE_class_union_bp']=='0' for r in feature_context if r['family']=='U2')
    assert gff_before_audit==(out/'Scotti_nuclear.ncRNA.combined.candidate.gff3').read_bytes()
    # Refuse overwriting a completed package; preserve original input bytes.
    before = trna.read_bytes(), protein.read_bytes()
    rerun = subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    assert rerun.returncode != 0 and 'Output exists' in rerun.stderr.decode()
    assert before==(trna.read_bytes(),protein.read_bytes())
print('PASS: selection/review separation, Rfam overlap marks, residual overlaps, E-value precision,')
print('      tRNA protection, minus-strand BED coordinates, union CDS overlap and overwrite protection.')
print('PASS: union repeat coverage, separation of small-RNA/TE/unknown classes, no candidate changes.')
