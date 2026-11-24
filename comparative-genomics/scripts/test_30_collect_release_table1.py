#!/usr/bin/env python3
"""Small synthetic regression tests; no real genomics data required."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).with_name("30_collect_release_table1.py")
spec = importlib.util.spec_from_file_location("collector", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class CollectorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def file(self, name, content):
        p = self.root/name
        p.write_text(content)
        return p

    def test_genome_stats(self):
        p = self.file("assembly.fa", ">a\n"+"A"*60+"\n>b\n"+"CG"*15+"\n>c\nNNNNNNNNNN\n")
        stats = m.fasta_stats(p, True)
        self.assertEqual((stats["genome_bp"], stats["n50_bp"], stats["n75_bp"], stats["l50"], stats["l75"]),
                         (100, 60, 30, 1, 2))
        self.assertAlmostEqual(stats["gc"], 100/3)
        self.assertEqual(stats["n_percent"], 10)
        q = self.file("wrapped.fa", ">c text\nnnnnn\nnnnnn\n>b\n"+"cg"*15+"\n>a\n"+"a"*60+"\n")
        self.assertEqual(stats["id_sequence_sha256"], m.fasta_stats(q, True)["id_sequence_sha256"])

    def test_duplicate_fasta(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            m.fasta_stats(self.file("duplicate.fa", ">a\nAA\n>a\nCC\n"))

    def test_annotation_and_representative_exons(self):
        gff = self.file("annotation.gff3", "##gff-version 3\n"
            "a\tt\tgene\t1\t100\t.\t+\t.\tID=gene1\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=rna-T1;Parent=gene1\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=e1;Parent=rna-T1\n"
            "a\tt\texon\t51\t100\t.\t+\t.\tID=e2;Parent=rna-T1\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=cds-P1;Parent=rna-T1;protein_id=P1\n")
        stats, problems = m.annotation_stats(gff, {"a": 100}, {"P1"})
        self.assertFalse(problems)
        self.assertEqual(stats["exons_per_gene"], 2)
        self.assertEqual(stats["mean_exon_bp"], 30)
        self.assertEqual(stats["mean_gene_bp"], 100)
        with self.assertRaisesRegex(ValueError, "coordinate mismatch"):
            m.annotation_stats(gff, {"a": 50}, {"P1"})

    def test_repeat_union(self):
        rm = self.file("assembly.fa.out", " 100 0.0 0.0 0.0 a 1 50 (50) + r DNA 1 50 (0) 1\n"
                       " 100 0.0 0.0 0.0 a 40 70 (30) + r DNA 1 31 (0) 2\n")
        self.assertEqual(m.repeat_percent(rm, {"a": 100}), 70)
        with self.assertRaisesRegex(ValueError, "length differs"):
            m.repeat_percent(rm, {"a": 101})

    def repeat_fixture(self):
        original = self.file("original.fa", ">a\n"+"A"*100+"\n>scaffold_61\nCCCCCCCCCC\n")
        selected = self.file("release.fa", ">a release\n"+"a"*50+"\n"+"a"*50+"\n")
        rm = self.file("old.fa.out", " 100 0.0 0.0 0.0 a 1 50 (50) + r DNA 1 50 (0) 1\n"
                       " 100 0.0 0.0 0.0 scaffold_61 1 10 (0) + r DNA 1 10 (0) 2\n")
        src = dict(genome=str(selected), repeatmasker_out=str(rm), repeatmasker_genome=str(original),
                   repeatmasker_exclude_ids=["scaffold_61"])
        return src, m.fasta_stats(selected, True), m.fasta_stats(original, True)

    def test_release_repeat_projection_exact_union_and_denominator(self):
        src, selected, original = self.repeat_fixture()
        result = m.repeat_evidence(src, selected, original)
        self.assertEqual(result["repeat_percent"], 50)
        self.assertEqual(result["masked_bp"], 50)
        self.assertEqual(result["excluded_masked_bp"], 10)
        self.assertEqual(result["excluded_genome_bp"], 10)
        self.assertEqual(result["denominator_bp"], 100)
        self.assertEqual(result["excluded_repeat_records"], 1)
        self.assertEqual(result["sequence_identity_status"], "all_retained_IDs_and_sequences_verified")

    def test_repeat_projection_rejects_changed_same_length_sequence(self):
        src, selected, original = self.repeat_fixture()
        changed = self.file("changed.fa", ">a\n"+"C"*100+"\n")
        with self.assertRaisesRegex(ValueError, "sequences changed"):
            m.repeat_evidence(src, m.fasta_stats(changed, True), original)

    def test_repeat_projection_requires_exact_exclusions_and_original(self):
        src, selected, original = self.repeat_fixture()
        with self.assertRaisesRegex(ValueError, "original input FASTA"):
            m.repeat_evidence(src, selected)
        src["repeatmasker_exclude_ids"] = []
        with self.assertRaisesRegex(ValueError, "ID difference"):
            m.repeat_evidence(src, selected, original)
        src["repeatmasker_exclude_ids"] = "scaffold_61"
        with self.assertRaisesRegex(ValueError, "must be a list"):
            m.repeat_evidence(src, selected, original)

    def test_excluded_repeat_coordinates_still_validated(self):
        src, selected, original = self.repeat_fixture()
        rm = Path(src["repeatmasker_out"])
        rm.write_text(rm.read_text().replace("scaffold_61 1 10 (0)", "scaffold_61 1 11 (0)"))
        with self.assertRaisesRegex(ValueError, "coordinates differ"):
            m.repeat_evidence(src, selected, original)

    def test_projection_allows_verified_zero_retained_repeat_hits(self):
        src, selected, original = self.repeat_fixture()
        Path(src["repeatmasker_out"]).write_text(
            " 100 0.0 0.0 0.0 scaffold_61 1 10 (0) + r DNA 1 10 (0) 2\n")
        self.assertEqual(m.repeat_evidence(src, selected, original)["repeat_percent"], 0)

    def test_actual_multilocularis_multipart_gene(self):
        spans = [(835414,835711), (835749,835921), (835959,836195),
                 (836232,836391), (841113,854231)]
        lines = [f"LN902841.1\tEMBL\tgene\t{a}\t{b}\t.\t+\t.\t"
                 "ID=gene-EmuJ_000747700.2;Name=EmuJ_000747700.2;gene_biotype=protein_coding\n"
                 for a,b in spans]
        lines += ["LN902841.1\tEMBL\tmRNA\t835414\t854231\t.\t+\t.\tID=T1;Parent=gene-EmuJ_000747700.2\n",
                  "LN902841.1\tEMBL\tmRNA\t841113\t854231\t.\t+\t.\tID=T1;Parent=gene-EmuJ_000747700.2\n",
                  "LN902841.1\tEMBL\tCDS\t835414\t835711\t.\t+\t0\tID=cds-P1;Parent=T1;protein_id=P1\n",
                  "LN902841.1\tEMBL\texon\t835414\t835711\t.\t+\t.\tID=E1;Parent=T1\n",
                  "LN902841.1\tEMBL\tgene\t900000\t900010\t.\t+\t.\tID=noncoding;gene_biotype=ncRNA\n",
                  "LN902841.1\tEMBL\tncRNA\t900000\t900010\t.\t+\t.\tID=R1;Parent=noncoding\n",
                  "LN902841.1\tEMBL\texon\t900000\t900010\t.\t+\t.\tID=RE1;Parent=R1\n"]
        p = self.file("multipart.gff3", "".join(lines))
        audit = []
        stats, problems = m.annotation_stats(p, {"LN902841.1": 1000000}, {"P1"}, audit)
        self.assertEqual(stats["genes"], 2)
        self.assertEqual(stats["coding_genes"], 1)
        self.assertEqual(stats["mean_gene_bp"], 18818)
        self.assertEqual(stats["coding_exons"], 1)
        self.assertEqual(stats["exons"], 2)
        self.assertEqual(stats["mean_exon_bp"], 298)
        self.assertFalse(problems)
        self.assertEqual(next(r for r in audit if r["kind"] == "gene" and r["rows"] > 1)["rows"], 5)

    def test_conflicting_duplicate_still_rejected(self):
        p = self.file("bad.gff3", "a\tt\tgene\t1\t10\t.\t+\t.\tID=g1\n"
                      "b\tt\tgene\t1\t10\t.\t+\t.\tID=g1\n")
        with self.assertRaisesRegex(ValueError, "Incompatible multipart"):
            m.annotation_stats(p, {"a": 100, "b": 100}, set())

    def test_gene_and_pseudogene_exon_mirrors_not_double_counted(self):
        # Small synthetic model of the observed D. latus hierarchy: coding and
        # pseudogene exon rows mirrored at mRNA level, plus one mRNA-only exon.
        lines = [
            "a\tt\tgene\t1\t100\t.\t+\t.\tID=g1\n",
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n",
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=t1;protein_id=P1\n",
            "a\tt\tgene\t201\t300\t.\t+\t.\tID=g2\n",
            "a\tt\tmRNA\t201\t300\t.\t+\t.\tID=t2;Parent=g2\n",
            "a\tt\tCDS\t201\t210\t.\t+\t0\tID=c2;Parent=t2;protein_id=P2\n",
            "a\tt\tpseudogene\t501\t600\t.\t+\t.\tID=pg\n",
            "a\tt\tmRNA\t501\t600\t.\t+\t.\tID=pt;Parent=pg\n",
            "a\tt\tCDS\t501\t510\t.\t+\t0\tID=pc;Parent=pt;protein_id=PP\n",
        ]
        for gene, transcript, intervals in [
                ("g1", "t1", [(1, 10), (51, 100)]),
                ("g2", "t2", [(201, 210), (251, 260)])]:
            for j, (start, end) in enumerate(intervals):
                for parent in (gene, transcript):
                    lines.append(f"a\tt\texon\t{start}\t{end}\t.\t+\t.\tID=e_{parent}_{j};Parent={parent}\n")
        lines += ["a\tt\texon\t501\t510\t.\t+\t.\tID=pe_gene;Parent=pg\n",
                  "a\tt\texon\t501\t510\t.\t+\t.\tID=pe_tx;Parent=pt\n",
                  "a\tt\texon\t551\t560\t.\t+\t.\tID=pe_only_tx;Parent=pt\n"]
        path = self.file("gene_mrna_mirrors.gff3", "".join(lines))
        before = path.read_bytes()
        audit = {}
        stats, problems = m.annotation_stats(path, {"a": 1000}, {"P1", "P2"}, exon_diagnostics=audit)
        self.assertFalse(problems)
        self.assertEqual((stats["genes"], stats["coding_genes"]), (2, 2))
        self.assertEqual((stats["exons"], stats["coding_exons"]), (6, 4))
        self.assertEqual((stats["mean_exon_bp"], stats["exons_per_gene"]), (20, 2))
        self.assertEqual(audit["raw_exon_records"], 11)
        self.assertEqual(audit["raw_gene_parent_links"], 4)
        self.assertEqual(audit["raw_pseudogene_parent_links"], 1)
        self.assertEqual(audit["raw_transcript_parent_links"], 6)
        self.assertEqual(audit["pseudogene_transcript_exon_associations"], 2)
        self.assertEqual(audit["pseudogene_transcripts"], 1)
        self.assertEqual(path.read_bytes(), before)

    def test_shared_and_repeated_exons_count_once_per_transcript(self):
        path = self.file("shared_exons.gff3", "a\tt\tgene\t1\t100\t.\t+\t.\tID=g\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=t1;protein_id=P1\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c2;Parent=t2;protein_id=P2\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=shared;Parent=t1,t2\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=repeated;Parent=t1\n"
            "a\tt\texon\t21\t30\t.\t+\t.\tID=extra;Parent=t1\n")
        audit = {}
        stats, problems = m.annotation_stats(path, {"a": 100}, {"P1"}, exon_diagnostics=audit)
        self.assertFalse(problems)
        self.assertEqual((stats["exons"], stats["coding_exons"]), (3, 3))
        self.assertEqual((stats["mean_exon_bp"], stats["exons_per_gene"]), (10, 2))
        self.assertEqual(audit["raw_exon_records"], 3)
        self.assertEqual(audit["raw_transcript_parent_links"], 4)
        self.assertEqual(audit["duplicate_transcript_parent_links_removed"], 1)

    def test_pseudogene_biotype_excluded_even_with_gene_feature(self):
        path = self.file("pseudo_biotype.gff3", "a\tt\tgene\t1\t100\t.\t+\t.\tID=g;gene_biotype=processed_pseudogene\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t;Parent=g\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c;Parent=t;protein_id=P\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=e;Parent=t\n")
        audit = {}
        stats, problems = m.annotation_stats(path, {"a": 100}, {"P"}, exon_diagnostics=audit)
        self.assertEqual(stats["genes"], 1)
        self.assertEqual(stats["coding_genes"], 0)
        self.assertEqual(stats["exons"], 1)
        self.assertIsNone(stats["coding_exons"])
        self.assertIsNone(stats["exons_per_gene"])
        self.assertTrue(problems)
        self.assertEqual(audit["pseudogene_transcripts"], 1)

    def test_missing_transcript_exons_not_replaced_by_gene_or_CDS(self):
        path = self.file("partial_exons.gff3", "a\tt\tgene\t1\t100\t.\t+\t.\tID=g\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=t1;protein_id=P1\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c2;Parent=t2;protein_id=P2\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=eg;Parent=g\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=e1;Parent=t1\n")
        audit = {}
        stats, _ = m.annotation_stats(path, {"a": 100}, {"P1"}, exon_diagnostics=audit)
        self.assertEqual(stats["exons"], 1)
        self.assertIsNone(stats["coding_exons"])
        self.assertIsNone(stats["mean_exon_bp"])
        self.assertEqual(audit["coding_transcripts_without_exons"], 1)
        self.assertFalse(audit["coding_exon_metrics_complete"])

    def test_CDS_only_annotation_does_not_invent_exons(self):
        path = self.file("cds_only.gff3", "a\tt\tgene\t1\t100\t.\t+\t.\tID=g\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c;Parent=g;protein_id=P\n")
        stats, problems = m.annotation_stats(path, {"a": 100}, {"P"})
        self.assertEqual(stats["coding_genes"], 1)
        for key in ("exons", "coding_exons", "mean_exon_bp", "exons_per_gene"):
            self.assertIsNone(stats[key])
        self.assertTrue(problems)

    def test_gene_level_exons_do_not_fill_missing_transcript_exons(self):
        path = self.file("gene_exons_only.gff3", "a\tt\tgene\t1\t100\t.\t+\t.\tID=g\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t;Parent=g\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c;Parent=t;protein_id=P\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=e;Parent=g\n")
        audit = {}
        stats, problems = m.annotation_stats(path, {"a": 100}, {"P"}, exon_diagnostics=audit)
        for key in ("exons", "coding_exons", "mean_exon_bp", "exons_per_gene"):
            self.assertIsNone(stats[key])
        self.assertEqual(audit["raw_gene_parent_links"], 1)
        self.assertEqual(audit["coding_transcripts_without_exons"], 1)
        self.assertTrue(problems)

    def test_idless_RNA_feature_does_not_invalidate_coding_metrics(self):
        path = self.file("anonymous_rna.gff3", "a\tt\tgene\t1\t100\t.\t+\t.\tID=g;pseudo=false\n"
            "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=t;Parent=g\n"
            "a\tt\tCDS\t1\t10\t.\t+\t0\tID=c;Parent=t;protein_id=P\n"
            "a\tt\texon\t1\t10\t.\t+\t.\tID=e;Parent=t\n"
            "a\tt\ttRNA\t200\t250\t.\t+\t.\tName=anonymous\n")
        stats, problems = m.annotation_stats(path, {"a": 1000}, {"P"})
        self.assertFalse(problems)
        self.assertEqual(stats["coding_genes"], 1)
        self.assertEqual(stats["coding_exons"], 1)

    def test_eleven_species_integration_and_wrong_lineage(self):
        genome = self.file("genome.fa", ">a\n"+"ACGT"*25+"\n")
        prot = self.file("protein.fa", ">g1.t1\nMPEPTIDE\n")
        gff = self.file("models.gff3", "a\tt\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                        "a\tt\tmRNA\t1\t100\t.\t+\t.\tID=g1.t1;Parent=g1\n"
                        "a\tt\tCDS\t1\t100\t.\t+\t0\tID=c1;Parent=g1.t1;protein_id=g1.t1\n"
                        "a\tt\texon\t1\t100\t.\t+\t.\tID=e_gene;Parent=g1\n"
                        "a\tt\texon\t1\t100\t.\t+\t.\tID=e1;Parent=g1.t1\n")
        manifest = m.default_sources(self.root)
        for row in manifest.values():
            row.update(genome=str(genome), primary=str(prot), isoforms=str(prot),
                       gff=str(gff), repeatmasker_out="")
        old = self.file("before_cleanup.fa", genome.read_text()+">scaffold_61\nAAAAAAAAAA\n")
        rm = self.file("before_cleanup.fa.out", " 100 0.0 0.0 0.0 a 1 50 (50) + r DNA 1 50 (0) 1\n"
                       " 100 0.0 0.0 0.0 scaffold_61 1 10 (0) + r DNA 1 10 (0) 2\n")
        manifest["S_solidus"].update(repeatmasker_out=str(rm), repeatmasker_genome=str(old))
        cfg = self.file("sources.json", json.dumps(manifest))
        audit = self.root/"audit"
        audit.mkdir()
        identities, reviews = [], []
        for sp in m.SPECIES:
            for mode, path, nuc in [("genome", genome, True), ("proteins", prot, False)]:
                identities.append(dict(species=sp, mode=mode, id_sequence_sha256=m.fasta_stats(path, nuc)["id_sequence_sha256"]))
                reviews.append(dict(species=sp, mode=mode, review_status="matched_current_input_check_run_time_provenance",
                                    C="80.0", S="79.0", D="1.0", F="5.0", M="15.0", n="672",
                                    dataset_date="2026-05-22", busco_version="6.0.0"))
        # One odb10 result must NOT populate the odb12 table.
        reviews[0]["n"] = "954"
        for name, rows in [("selected_inputs.tsv", identities),
                           ("eleven_species.genome_and_protein.review.tsv", reviews)]:
            with (audit/name).open("w", newline="") as f:
                w = csv.DictWriter(f, rows[0].keys(), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
        out = self.root/"results"
        subprocess.run([sys.executable, str(SCRIPT), "--root", str(self.root), "--inputs", str(cfg),
                        "--audit-dir", str(audit), "--out", str(out)], check=True, capture_output=True)
        table = {r["Assembly / annotation feature"]: r for r in m.read_tsv(out/"Table1_11species.tsv")}
        self.assertEqual(table["Complete genome BUSCO (%)"]["S_solidus"], "n.a.")
        self.assertEqual(table["Complete protein BUSCO (%)"]["S_solidus"], "80.0")
        self.assertEqual(table["Mean exons per gene (representative isoform)"]["S_cotti"], "1.00")
        self.assertEqual(table["Number of transposon proteins (representative set)"]["S_cotti"], "n.a.")
        self.assertEqual(table["Bases masked by RepeatMasker (%)"]["S_solidus"], "50.00")
        repeat_audit = m.read_tsv(out/"repeatmasker_release_audit.tsv")
        self.assertEqual(repeat_audit[0]["excluded_ids"], "scaffold_61")
        self.assertEqual(repeat_audit[0]["masked_bp"], "50")
        exon_audit = m.read_tsv(out/"exon_transcript_audit.tsv")
        self.assertEqual(len(exon_audit), 11)
        self.assertEqual(exon_audit[0]["method"], m.EXON_METHOD)
        self.assertEqual(exon_audit[0]["coding_transcript_exon_associations_observed"], "1")
        self.assertEqual(exon_audit[0]["raw_exon_records"], "2")
        self.assertEqual(exon_audit[0]["raw_gene_parent_links"], "1")
        self.assertEqual(table["Exon-transcript associations (protein-coding transcripts)"]["S_cotti"], "1")
        self.assertEqual(json.loads((out/"collection.json").read_text())["exon_method"], m.EXON_METHOD)
        self.assertIn(">scaffold_61", old.read_text())
        self.assertIn("scaffold_61", rm.read_text())
        self.assertEqual(len(next(iter(table.values()))), 12)
        # Existing output directories are never overwritten.
        result = subprocess.run([sys.executable, str(SCRIPT), "--root", str(self.root), "--inputs", str(cfg),
                                 "--out", str(out)], capture_output=True)
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
