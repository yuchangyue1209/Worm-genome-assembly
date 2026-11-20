"""Small synthetic checks for overlap arithmetic and non-destructive triage."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv, importlib.util, json, subprocess, sys

MODULE = Path(__file__).with_name('32_build_scotti_te_review.py')
spec = importlib.util.spec_from_file_location('build_scotti_te_review', MODULE)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)
Coverage, classify = module.Coverage, module.classify

idx=Coverage([(0,100),(50,125),(225,275)])
assert idx.overlap(0,150)+idx.overlap(200,350)==175
assert idx.overlap(125,225)==0
assert idx.overlap(100,250)==50
assert classify({'host_context':{'TERT'},'RT_or_integrase':{'RT'}},100,0)=='host_context_review'
assert classify({},100,0)=='repeat_context_only'
assert classify({},0,0)=='no_screening_signal'

with TemporaryDirectory() as temp:
    root=Path(temp);c=root/'candidate';e=c/'functional_evidence_v1';e.mkdir(parents=True)
    ids=[f'g{i}.t1' for i in range(1,8)]
    for name in ['Scotti_nuclear.filtered7.proteins.fa','Scotti_nuclear.filtered7.longest.proteins.fa']:
        (c/name).write_text(''.join(f'>{p}\n'+('A'*100)+'\n' for p in ids))
    rows=[]
    for i,p in enumerate(ids,1):
        rows.extend([f'c{i}\tfixture\tgene\t1\t350\t.\t+\t.\tID=g{i}',
                     f'c{i}\tfixture\tmRNA\t1\t350\t.\t+\t.\tID={p};Parent=g{i}',
                     f'c{i}\tfixture\tCDS\t1\t150\t.\t+\t0\tParent={p}',
                     f'c{i}\tfixture\tCDS\t201\t350\t.\t+\t0\tParent={p}'])
    (c/'Scotti_nuclear.filtered7.candidate.gff3').write_text('##gff-version 3\n'+'\n'.join(rows)+'\n')
    fields=['protein_id','gene_id','length_aa','swissprot_title','swissprot_evalue','swissprot_query_coverage_pct','eggnog_description','interpro_descriptions','eggnog_PFAMs']
    with (e/'Scotti_filtered7.all_isoforms.evidence.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for i,p in enumerate(ids,1):
            r=dict.fromkeys(fields,'');r.update(protein_id=p,gene_id=f'g{i}',length_aa=100)
            if i==1:r['swissprot_title']='Retrotransposon Pol polyprotein'
            if i==2:r['swissprot_title']='Telomerase reverse transcriptase'
            w.writerow(r)
    (e/'Scotti_filtered7.interproscan.tsv').write_text('\n'.join([
        'g1.t1\tmd5\t100\tPfam\tPF00000\tReverse transcriptase\t10\t90\t1e-20\tT\tdate\t-\t-',
        'g3.t1\tmd5\t100\tPfam\tPF00001\tRibonuclease H\t10\t90\t1e-10\tT\tdate\t-\t-'])+'\n')
    rm=root/'repeat.out'
    lines=[]
    for seq,start,end,cls in [('c1',1,100,'LINE/CR1'),('c1',51,125,'LINE/CR1'),('c1',226,275,'Unknown'),('c1',301,350,'tRNA'),('c4',1,350,'DNA/TcMar'),('c5',1,350,'Unknown'),('c6',1,350,'rRNA')]:
        lines.append(f'100 1.0 0.0 0.0 {seq} {start} {end} (0) + repeat {cls} 1 100 (0) 1')
    rm.write_text('\n'.join(lines)+'\n')
    out=root/'review'
    subprocess.run([sys.executable,str(MODULE),
        '--candidate',str(c),'--repeat-out',str(rm),'--out',str(out),
        '--expected-proteins','7','--expected-genes','7'],check=True,stdout=subprocess.PIPE,text=True)
    with (out/'Scotti_TE_review.all_isoforms.tsv').open() as f:
        results={r['protein_id']:r for r in csv.DictReader(f,delimiter='\t')}
    assert float(results['g1.t1']['CDS_RM_all_pct'])==75
    assert abs(float(results['g1.t1']['CDS_RM_classified_TE_pct'])-41.6667)<.0001
    assert results['g1.t1']['review_tier']=='priority_functional_review'
    assert results['g2.t1']['review_tier']=='host_context_review'
    assert results['g3.t1']['review_tier']=='generic_domain_review'
    assert results['g4.t1']['review_tier']=='repeat_context_only'
    assert results['g5.t1']['review_tier']=='repeat_context_only'
    assert results['g6.t1']['review_tier']=='no_screening_signal'
    assert results['g7.t1']['review_tier']=='no_screening_signal'
    assert (out/'Scotti.primary.TE_functional_candidate_ids.txt').read_text().splitlines()==['g1.t1','g2.t1']
    assert not list(out.glob('*.fa'))
    assert json.loads((out/'summary.json').read_text())['genes']==7
    # Original protein inputs remain byte-for-byte unchanged.
    assert len(list(out.glob('*review*.tsv')))>=4
print('PASS: interval union/CDS boundaries, host-context flag, functional/domain/repeat-only tiers, candidate IDs, and no sequence filtering.')
