"""Exercise all ten missing jobs and a verified resume using a simulated BUSCO."""
from pathlib import Path
import csv
import json
import sys
import tempfile
from types import SimpleNamespace
import audit_final_busco_inputs as audit
import run_missing_genome_busco as runner

with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    project = root/'project'
    input_dir = root/'inputs'
    input_dir.mkdir()
    audit_dir = project/'03_qc/final_busco_audit_20261120_123456'
    audit_dir.mkdir(parents=True)
    inputs, review, reports = [], [], []
    for species in audit.SPECIES:
        for mode in ['genome','proteins']:
            fasta = input_dir/f'{species}.{mode}.fa'
            fasta.write_text('>id1\nACGTA\n>id2\nTTT\n')
            row = dict(species=species,mode=mode,selected_input=str(fasta),input_status='present',
                       **audit.fingerprint(fasta,{}))
            inputs.append(row)
            matched = mode=='proteins' or species=='S_solidus'
            review.append({**row,'review_status':runner.MATCHED if matched else 'NO_MATCHING_ODB12_REPORT',
                           'C':'73.5' if matched else '', 'S':'72.5' if matched else '',
                           'D':'1.0' if matched else '', 'F':'9.8' if matched else '',
                           'M':'16.7' if matched else '', 'n':'672' if matched else '',
                           'busco_version':'6.0.0' if matched else '',
                           'dataset_date':'2026-05-22' if matched else ''})
            if matched:
                reports.append(dict(species=species,mode=mode,pipeline='euk_genome_min' if mode=='genome' else 'proteins',
                                    lineage='metazoa_odb12',match_status='current_ID_sequences_match_run_time_hash_not_proven',
                                    summary=str(root/'no_such_summary.txt')))
    audit.tsv(audit_dir/'selected_inputs.tsv', inputs, list(inputs[0]))
    audit.tsv(audit_dir/'eleven_species.genome_and_protein.review.tsv',review,list(review[0]))
    audit.tsv(audit_dir/'all_existing_busco_reports.tsv',reports,list(reports[0]))
    db = root/'metazoa_odb12'
    db.mkdir()
    (db/'dataset.cfg').write_text('creation_date=2026-05-22\nnumber_of_BUSCOs=672\n')
    (db/'hmms').mkdir()
    (db/'hmms/fake.hmm').write_text('fake fixture\n')
    (db/'scores_cutoff').write_text('fixture\n')
    (db/'refseq_db.faa.gz').write_bytes(b'fixture')
    # Regression: original odb12 deliberately lacks lengths_cutoff.
    assert not (db/'lengths_cutoff').exists()
    assert runner.dataset_ok(db), runner.dataset_problems(db)
    assert runner.find_database(db, project, audit_dir) == db.resolve()
    score_file = db/'scores_cutoff'
    score_file.rename(db/'scores_cutoff.saved')
    assert not runner.dataset_ok(db)
    assert 'missing/empty scores_cutoff' in runner.dataset_problems(db)
    try:
        runner.find_database(db, project, audit_dir)
    except ValueError as error:
        assert 'missing/empty scores_cutoff' in str(error)
    else:
        raise AssertionError('A missing required scores_cutoff was not rejected')
    (db/'scores_cutoff.saved').rename(score_file)
    assert runner.dataset_ok(db)
    calls = []
    def simulated(command, **kwargs):
        if '--version' in command:
            return SimpleNamespace(stdout='BUSCO 6.0.0\n',stderr='')
        if '--help' in command:
            return SimpleNamespace(stdout='--miniprot --offline\n',stderr='')
        name = command[command.index('-o')+1]
        fasta = command[command.index('-i')+1]
        folder = Path(command[command.index('--out_path')+1])/name
        folder.mkdir(parents=True)
        summary = folder/f'short_summary.specific.metazoa_odb12.{name}.txt'
        summary.write_text('# BUSCO version is: 6.0.0\n'
                           '# The lineage dataset is: metazoa_odb12 (Creation date: 2026-05-22)\n'
                           f'# Summarized benchmarking in BUSCO notation for file {fasta}\n'
                           '# BUSCO was run in mode: euk_genome_min\n'
                           'C:73.5%[S:72.5%,D:1.0%],F:9.8%,M:16.7%,n:672\n')
        calls.append(name)
        return SimpleNamespace(returncode=0)
    runner.subprocess.run = simulated
    runner.shutil.which = lambda name: '/fixture/bin/busco'
    out = root/'run'
    sys.argv = [str(Path(runner.__file__)),'--root',str(project),'--audit-dir',str(audit_dir),
                '--lineage-db',str(db),'--out',str(out),'--cpu','2']
    runner.main()
    assert len(calls)==10 and 'S_solidus' not in calls
    assert len(list((out/'inputs').glob('*.fa')))==10
    with (out/'eleven_species.genome_and_protein.odb12.tsv').open() as handle:
        rows = list(csv.DictReader(handle,delimiter='\t'))
    assert len(rows)==22
    assert sum(r['input_provenance']=='frozen_input_run_hash_verified' for r in rows)==10
    sys.argv.append('--resume')
    runner.main()
    assert len(calls)==10
    # An audited input mutation must block a subsequent resume before any job.
    Path(inputs[0]['selected_input']).write_text('>id1\nAAAAA\n')
    try:
        runner.main()
    except ValueError as error:
        assert 'Input changed after audit' in str(error)
    else:
        raise AssertionError('Changed input was not rejected')
    assert len(calls)==10
print('PASS: odb12 without lengths_cutoff, missing required file diagnosis, ten simulated genome jobs, frozen provenance, 22-row comparison, verified resume, changed-input rejection.')
