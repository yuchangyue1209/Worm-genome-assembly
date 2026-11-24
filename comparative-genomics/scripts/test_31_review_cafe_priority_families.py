#!/usr/bin/env python3
"""Synthetic regression tests; no real server files or external tools needed."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).with_name("31_review_cafe_priority_families.py")
spec = importlib.util.spec_from_file_location("review", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def file(self, name, text):
        path = self.root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def hit(self, pid, seq, signature, start, end, ipr="-", analysis="Pfam"):
        return "\t".join([pid, hashlib.md5(seq.encode()).hexdigest(), str(len(seq)), analysis,
            signature, "Test feature", str(start), str(end), "1e-20", "T", "2026-11-23", ipr, "Test InterPro"])+"\n"

    def fixture(self):
        ogs = m.DEFAULT_OGS
        table = "Orthogroup\tS_cotti\tS_solidus\n"+"".join(
            f"{og}\tg{i}.t1\tg{i}.t1\n" for i, og in enumerate(ogs, 1))
        self.results = self.root/"results"
        self.file("results/Orthogroups/Orthogroups.tsv", table)
        self.seqs = {f"g{i}.t1": "M"+"A"*(119+i*10) for i in range(1, 5)}
        for sp in m.FOCAL:
            self.file(f"input/{sp}.faa", "".join(f">{pid}\n{seq}\n" for pid, seq in self.seqs.items())+
                      ">outside_same_sequence\n"+self.seqs["g1.t1"]+"\n")
            gff = []
            ipr = []
            for i, (pid, seq) in enumerate(self.seqs.items(), 1):
                start, end = i*2000+1, i*2000+3*len(seq)
                tid, gid = f"rna-model{i}", f"gene{i}"
                gff += [f"chr1\tt\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}\n",
                        f"chr1\tt\tmRNA\t{start}\t{end}\t.\t+\t.\tID={tid};Parent={gid};protein_id={pid}\n",
                        f"chr1\tt\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds-{pid};Parent={tid};protein_id={pid}\n"]
                if i == 1:
                    ipr.append(self.hit(pid, seq, "PF03645", 1, 100, "IPR005334"))
                elif i == 2:
                    ipr.append(self.hit(pid, seq, "PF04968", 1, 50, "IPR007051"))
                    if sp == "S_cotti":
                        ipr.append(self.hit(pid, seq, "PF04969", 40, 100, "IPR007052"))
                elif i == 3:
                    ipr.append(self.hit(pid, seq, "PF00268", 1, 100, "IPR033909"))
                    if sp == "S_solidus":
                        ipr.append(self.hit(pid, seq, "PS_TEST", 30, 30, "IPR030475", "ProSitePatterns"))
                elif sp == "S_solidus":
                    ipr.append(self.hit(pid, seq, "PF00650", 1, 100, "IPR001251"))
            self.file(f"{sp}.gff3", "".join(gff))
            self.file(f"{sp}.tsv", "".join(ipr))
        self.out = self.root/"output"
        self.command = [sys.executable, str(SCRIPT), "--results", str(self.results),
            "--input-dir", str(self.root/"input"), "--out", str(self.out), "--expected-species", "2"]
        for sp in m.FOCAL:
            self.command += ["--gff", f"{sp}={self.root/(sp+'.gff3')}",
                             "--interpro", f"{sp}={self.root/(sp+'.tsv')}"]

    def test_union_counts_overlap_once(self):
        self.assertEqual(m.union_bp([(1, 50), (40, 100), (1, 50), (101, 110)]), 110)

    def test_duplicate_fasta_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            m.read_fasta(self.file("bad.fa", ">x\nMA\n>x\nMC\n"))

    def test_no_gene_id_guessing(self):
        path = self.file("absent.gff3", "chr1\tt\tgene\t1\t30\t.\t+\t.\tID=g1\n")
        loc, cds = m.locations(path, {"g1.t1"})
        self.assertEqual(loc["g1.t1"]["Mapping_status"], "Missing_or_ambiguous_transcript")
        self.assertNotIn("Gene_ID", loc["g1.t1"])
        self.assertFalse(cds)

    def test_explicit_mapping_and_bed_coordinates(self):
        self.fixture()
        loc, cds = m.locations(self.root/"S_cotti.gff3", set(self.seqs))
        self.assertEqual(loc["g1.t1"]["Gene_ID"], "gene1")
        self.assertEqual(loc["g1.t1"]["Start"], 2001)
        self.assertEqual(loc["g1.t1"]["CDS_union_bp"], len(self.seqs["g1.t1"])*3)
        self.assertEqual(cds["g1.t1"][0][1], 2001)

    def test_md5_mismatch_stops_before_output(self):
        self.fixture()
        path = self.root/"S_cotti.tsv"
        path.write_text(path.read_text().replace(hashlib.md5(self.seqs["g1.t1"].encode()).hexdigest(), "0"*32))
        result = subprocess.run(self.command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("sequence mismatch", result.stderr)
        self.assertFalse(self.out.exists())

    def test_feature_outside_protein_rejected(self):
        seq = "MPEPTIDE"
        path = self.file("outside.tsv", self.hit("p", seq, "PF_TEST", 1, 99))
        with self.assertRaisesRegex(ValueError, "outside protein"):
            m.annotations(path, {"p": seq})

    def test_single_terminal_stop_is_verified_and_recorded(self):
        seq = "M"+"A"*79
        path = self.file("terminal.tsv", self.hit("p", seq, "PF03645", 1, 80))
        hits = m.annotations(path, {"p": seq+"*"})
        self.assertEqual(hits["p"][0]["Sequence_match_rule"], "Single_terminal_stop_removed")
        self.assertEqual(hits["p"][0]["InterPro_sequence_length"], 80)
        exact = m.annotations(path, {"p": seq})
        self.assertEqual(exact["p"][0]["Sequence_match_rule"], "Exact_input_sequence")

    def test_terminal_stop_does_not_hide_internal_stops_or_substitutions(self):
        seq = "MPEPTIDE"
        path = self.file("strict.tsv", self.hit("p", seq, "PF_TEST", 1, 8))
        for changed in ["MPEP*TIDE*", seq+"**", "MPEPTIDF*"]:
            with self.subTest(sequence=changed):
                with self.assertRaisesRegex(ValueError, "sequence mismatch"):
                    m.annotations(path, {"p": changed})

    def test_terminal_stop_length_is_also_checked(self):
        seq = "MPEPTIDE"
        path = self.file("wrong_length.tsv", self.hit("p", seq, "PF_TEST", 1, 8))
        path.write_text(path.read_text().replace("\t8\tPfam", "\t9\tPfam"))
        with self.assertRaisesRegex(ValueError, "sequence mismatch"):
            m.annotations(path, {"p": seq+"*"})

    def test_terminal_stop_end_to_end_preserves_fasta_and_uses_aa_length(self):
        self.fixture()
        path = self.root/"input/S_cotti.faa"
        path.write_text(path.read_text().replace(">g1.t1\n"+self.seqs["g1.t1"]+"\n",
            ">g1.t1\n"+self.seqs["g1.t1"]+"*\n"))
        before = m.file_sha256(path)
        result = subprocess.run(self.command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        _, rows = m.read_table(self.out/"priority_protein_review.tsv")
        row = next(r for r in rows if (r["Species"], r["Protein_ID"]) == ("S_cotti", "g1.t1"))
        length = len(self.seqs["g1.t1"])
        self.assertEqual(row["Length_aa"], str(length))
        self.assertEqual(row["Input_length_including_terminal_stop"], str(length+1))
        self.assertEqual(row["Terminal_stop_symbols"], "1")
        self.assertEqual(row["InterPro_sequence_match"], "Single_terminal_stop_removed")
        self.assertEqual(float(row["Pfam_coverage_percent"]), round(10000/length, 2))
        self.assertNotIn("CDS_length_review_phase_partiality", row["Flags"])
        audit = m.read_table(self.out/"protein_sequence_match_audit.tsv")[1]
        audited = next(r for r in audit if (r["Species"], r["Protein_ID"]) == ("S_cotti", "g1.t1"))
        self.assertNotEqual(audited["Input_sequence_MD5"], audited["Amino_acid_sequence_MD5"])
        self.assertEqual(before, m.file_sha256(path))
        exported = m.read_fasta(self.out/"all_species_fastas/OG0000875.faa")
        self.assertTrue(exported["S_cotti|g1.t1"].endswith("*"))

    def test_conflicting_cds_coordinates_rejected(self):
        self.fixture()
        path = self.root/"S_cotti.gff3"
        path.write_text(path.read_text()+"chr2\tt\tCDS\t2001\t2010\t.\t+\t0\tParent=rna-model1\n")
        with self.assertRaisesRegex(ValueError, "Conflicting gene/transcript/CDS"):
            m.locations(path, set(self.seqs))

    def test_integration_and_source_preservation(self):
        self.fixture()
        source_hashes = {path: m.file_sha256(path) for path in self.root.rglob("*") if path.is_file()}
        result = subprocess.run(self.command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Focal proteins reviewed: 8", result.stdout)
        _, rows = m.read_table(self.out/"priority_protein_review.tsv")
        keyed = {(r["Species"], r["Protein_ID"]): r for r in rows}
        cotti = keyed["S_cotti", "g2.t1"]
        self.assertEqual(cotti["All_expected_Pfam_detected"], "1")
        self.assertEqual(cotti["Pfam_covered_aa"], "100")
        self.assertEqual(keyed["S_solidus", "g2.t1"]["All_expected_Pfam_detected"], "0")
        self.assertEqual(keyed["S_cotti", "g3.t1"]["RNR_active_site_entry_detected"], "0")
        self.assertEqual(keyed["S_solidus", "g3.t1"]["RNR_active_site_entry_detected"], "1")
        missing = keyed["S_cotti", "g4.t1"]
        self.assertEqual(missing["Annotation_status"], "No_record_in_InterPro_TSV")
        self.assertEqual(missing["Pfam_coverage_percent"], "NA")
        self.assertEqual(missing["All_expected_Pfam_detected"], "NA")
        self.assertEqual(keyed["S_cotti", "g1.t1"]["Identical_sequence_count_in_species_input"], "2")
        self.assertIn("chr1\t2000\t", (self.out/"S_cotti.priority_cds.bed").read_text())
        self.assertEqual(len(m.read_fasta(self.out/"all_species_fastas/OG0000875.faa")), 2)
        self.assertEqual(len(m.read_table(self.out/"family_copy_counts_allspecies.tsv")[1]), 4)
        provenance = json.loads((self.out/"provenance.json").read_text())
        self.assertFalse(provenance["RNA_support_tested"])
        for path, digest in source_hashes.items():
            self.assertEqual(m.file_sha256(path), digest)
        again = subprocess.run(self.command, capture_output=True, text=True)
        self.assertNotEqual(again.returncode, 0)
        self.assertIn("already exists", again.stderr)


if __name__ == "__main__":
    unittest.main()
