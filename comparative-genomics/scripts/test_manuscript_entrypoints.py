#!/usr/bin/env python3
"""Synthetic checks only; no biological jobs or release files are touched."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parent
REPO = SCRIPTS.parents[1]
spec = importlib.util.spec_from_file_location('verified_ipr', SCRIPTS/'interpro_verified_run.py')
ipr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ipr)


class ManuscriptSafetyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bin = self.root/'bin'
        self.bin.mkdir()
        self.env = os.environ.copy()
        self.env['PATH'] = str(self.bin)+os.pathsep+str(Path(sys.executable).parent)+os.pathsep+self.env['PATH']
        self.make_exec('sha256sum', 'import hashlib,sys\nfrom pathlib import Path\nfor a in sys.argv[1:]: print(hashlib.sha256(Path(a).read_bytes()).hexdigest(), a)\n')

    def make_exec(self, name, code):
        path = self.bin/name
        path.write_text('#!'+sys.executable+'\n'+code)
        path.chmod(0o755)
        return path

    def run_shell(self, name):
        return subprocess.run(['bash', str(SCRIPTS/name)], env=self.env, capture_output=True, text=True)

    def staging_fixture(self):
        for name, data in [('genome.fa', '>chr1\nACGT\n'), ('genes.gff', '##gff-version 3\n'),
                           ('primary.fa', '>g1.t1\nMA*\n'), ('all.fa', '>g1.t1\nMA*\n>g1.t2\nMAC*\n')]:
            (self.root/name).write_text(data)
        self.env.update(STD=str(self.root/'standard'), COTTI_GENOME=str(self.root/'genome.fa'),
                        COTTI_GFF=str(self.root/'genes.gff'), COTTI_PROTEINS=str(self.root/'primary.fa'),
                        COTTI_ALL_PROTEINS=str(self.root/'all.fa'), COTTI_CDS='',
                        EXPECTED_PROTEINS='1', EXPECTED_ALL_PROTEINS='2',
                        EXPECTED_SHA256=hashlib.sha256((self.root/'primary.fa').read_bytes()).hexdigest())

    def test_stage_distinct_raw_primary_and_idempotence(self):
        self.staging_fixture()
        for _ in range(2):
            r = self.run_shell('00d_stage_scotti_candidate.sh')
            self.assertEqual(r.returncode, 0, r.stderr)
        std = self.root/'standard'
        self.assertEqual((std/'proteins_primary/S_cotti.faa').read_bytes(), (self.root/'primary.fa').read_bytes())
        self.assertEqual((std/'proteins_raw/S_cotti.faa').read_bytes(), (self.root/'all.fa').read_bytes())

    def test_conflicting_raw_link_stops_before_other_links(self):
        self.staging_fixture()
        raw = self.root/'standard/proteins_raw/S_cotti.faa'
        raw.parent.mkdir(parents=True)
        raw.symlink_to(self.root/'primary.fa')
        before = os.readlink(raw)
        r = self.run_shell('00d_stage_scotti_candidate.sh')
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('Conflicting destination', r.stderr)
        self.assertEqual(os.readlink(raw), before)
        self.assertFalse((self.root/'standard/genomes').exists())

    def test_missing_optional_cds_stops_before_staging(self):
        self.staging_fixture()
        self.env['COTTI_CDS'] = str(self.root/'missing.fa')
        self.assertNotEqual(self.run_shell('00d_stage_scotti_candidate.sh').returncode, 0)
        self.assertFalse((self.root/'standard').exists())

    def test_primary_sequence_must_match_raw(self):
        self.staging_fixture()
        (self.root/'all.fa').write_text('>g1.t1\nMC*\n>g1.t2\nMAC*\n')
        self.assertNotEqual(self.run_shell('00d_stage_scotti_candidate.sh').returncode, 0)
        self.assertFalse((self.root/'standard').exists())

    def test_raw_duplicate_ids_rejected(self):
        self.staging_fixture()
        (self.root/'all.fa').write_text('>g1.t1\nMA*\n>g1.t1\nMAC*\n')
        self.assertNotEqual(self.run_shell('00d_stage_scotti_candidate.sh').returncode, 0)
        self.assertFalse((self.root/'standard').exists())

    def busco_fixture(self):
        for lineage, n in [('metazoa_odb12',672), ('metazoa_odb10',954)]:
            path = self.root/lineage/'S_cotti'/'short_summary.txt'
            path.parent.mkdir(parents=True)
            path.write_text(f'C:73.5%[S:72.5%,D:1.0%],F:9.8%,M:16.7%,n:{n}\n')

    def test_busco_defaults_include_cotti_and_only_odb12(self):
        self.busco_fixture()
        out = self.root/'summary.tsv'
        r = subprocess.run([sys.executable,str(SCRIPTS/'02e_collect_busco_summaries.py'),
                            '--root',str(self.root),'--out',str(out)], capture_output=True,text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn('S_cotti\tmetazoa_odb12', out.read_text())
        self.assertNotIn('metazoa_odb10', out.read_text())

    def test_busco_historical_lineages_explicit(self):
        self.busco_fixture()
        out = self.root/'summary.tsv'
        subprocess.run([sys.executable,str(SCRIPTS/'02e_collect_busco_summaries.py'),
                        '--root',str(self.root),'--out',str(out),'--lineages',
                        'metazoa_odb10','metazoa_odb12'],check=True,capture_output=True)
        self.assertEqual(len(out.read_text().splitlines()),3)

    def ipr_fixture(self, fail=False, empty=False):
        protein = self.root/'input.fa'
        protein.write_text('>g1.t1 description\nMA*\n>g2.t1\nMC\n')
        exe = self.make_exec('interproscan',
            'import sys,hashlib\nfrom pathlib import Path\n'
            'a=sys.argv; p=Path(a[a.index("-i")+1]); o=Path(a[a.index("-b")+1]+".tsv")\n'
            's=p.read_text().splitlines(); rows=[]\n'
            'for i in range(0,len(s),2):\n'
            ' rows.append("\\t".join([s[i][1:],hashlib.md5(s[i+1].encode()).hexdigest(),str(len(s[i+1])),"Pfam","PF1","test","1","2","1e-5","T","today"]))\n'
            + ('o.write_text("")\n' if empty else 'o.write_text("\\n".join(rows)+"\\n")\n')
            + ('sys.exit(1)\n' if fail else ''))
        return protein, self.root/'out/S_cotti', exe, 1, self.root/'logs/S_cotti.log'

    def test_interpro_normalizes_only_terminal_stop(self):
        args = self.ipr_fixture()
        original = args[0].read_bytes()
        ipr.run(*args)
        self.assertEqual(args[0].read_bytes(),original)
        self.assertNotIn('*',Path(str(args[1])+'.search.faa').read_text())
        audit=json.loads(Path(str(args[1])+'.input.json').read_text())
        self.assertTrue(audit['proteins'][0]['terminal_stop_removed'])
        ipr.run(*args)  # verified resume, not a second analysis

    def test_interpro_legacy_tsv_not_overwritten_or_rerun(self):
        args=self.ipr_fixture()
        path=Path(str(args[1])+'.tsv'); path.parent.mkdir(); path.write_text('old evidence\n')
        with self.assertRaises(FileExistsError): ipr.run(*args)
        self.assertEqual(path.read_text(),'old evidence\n')
        self.assertFalse(args[-1].exists())

    def test_interpro_failure_has_no_completion_marker(self):
        args=self.ipr_fixture(fail=True)
        with self.assertRaises(subprocess.CalledProcessError): ipr.run(*args)
        self.assertFalse(Path(str(args[1])+'.complete.json').exists())

    def test_interpro_tampered_output_rejected(self):
        args=self.ipr_fixture(); ipr.run(*args)
        Path(str(args[1])+'.tsv').write_text('changed\n')
        with self.assertRaises(ValueError): ipr.run(*args)

    def test_interpro_empty_successful_output_is_allowed(self):
        args=self.ipr_fixture(empty=True); ipr.run(*args)
        self.assertEqual(json.loads(Path(str(args[1])+'.complete.json').read_text())['proteins_with_tsv_rows'],0)

    def test_interpro_invalid_sequences_fail_before_output(self):
        for sequence in ['MA*X','MA**','', 'ma']:
            args=self.ipr_fixture(); args[0].write_text('>g1\n'+sequence+'\n')
            with self.assertRaises(ValueError): ipr.run(*args)
            self.assertFalse(args[1].parent.exists())

    def test_interpro_row_md5_checked(self):
        path=self.root/'row.tsv'
        path.write_text('g1\tbad\t2\tPfam\tPF1\ttest\t1\t2\t0\tT\ttoday\n')
        with self.assertRaises(ValueError): ipr.validate_tsv(path, {'g1':'MA'})

    def test_interpro_shell_uses_derived_input(self):
        protein,prefix,exe,threads,log=self.ipr_fixture()
        std=self.root/'std/proteins_primary'; std.mkdir(parents=True)
        (std/'S_cotti.faa').write_bytes(protein.read_bytes())
        config=self.root/'config.sh'
        config.write_text(f'STD="{std.parent}"\nFUNC="{self.root}/functional"\nINTERPROSCAN="{exe}"\nTHREADS=1\nSPECIES=(S_cotti)\n')
        self.env['CONFIG']=str(config)
        r=self.run_shell('04_interproscan.sh')
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertTrue((self.root/'functional/interpro/S_cotti.complete.json').is_file())

    def test_tagseq_empty_and_nonempty_exclusion_lists(self):
        script=(REPO/'solidus/workflow-scripts/21_tagseq_3prime_refinement.sh').read_text()
        expression=re.search(r"awk -F'\\t' '([^']*bad\[\$1\][^']*)'",script).group(1)
        ids=self.root/'bad.ids'; bed=self.root/'input.bed'
        bed.write_text('chr1\t0\t1\tkeep\nchr1\t2\t3\tremove\n')
        for excluded, expected in [('',2),('remove\n',1)]:
            ids.write_text(excluded)
            r=subprocess.run(['awk','-F','\t',expression,str(ids),str(bed)],capture_output=True,text=True,check=True)
            self.assertEqual(len(r.stdout.splitlines()),expected)

    def test_tagseq_missing_gnu_awk_stops_before_jobs(self):
        for name in ['cutadapt','samtools','bedtools','sort','gzip','zcat','md5sum','awk']:
            self.make_exec(name,'import sys\nsys.exit(0)\n')
        script=REPO/'solidus/workflow-scripts/21_tagseq_3prime_refinement.sh'
        r=subprocess.run(['bash',str(script)],env=self.env,capture_output=True,text=True)
        self.assertNotEqual(r.returncode,0)
        self.assertIn('GNU awk',r.stderr)


if __name__=='__main__':
    unittest.main()
