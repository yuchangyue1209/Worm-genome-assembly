#!/usr/bin/env python3
"""Synthetic regression fixtures, not a replacement for real Figure 1 data."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

SCRIPT = Path(__file__).with_name("17_plot_orthology_summary_clean.py")
spec = importlib.util.spec_from_file_location("orthology_plot", SCRIPT)
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)
TREE = "(D_latus:0.05,(((((E_granulosus:0.01,E_multilocularis:0.02)100/100:0.06,(T_multiceps:0.025,T_solium:0.019)100/100:0.06)100/100:0.14,H_microstoma:0.35)100/100:0.70,(S_proliferum:0.025,Spirometra_Aus1:0.02)100/100:0.05)100/100:0.019,(S_cotti:0.01,S_solidus:0.017)100/100:0.04)100/100:0.10,L_intestinalis:0.054);"
CYCLO = {"E_granulosus", "E_multilocularis", "T_multiceps", "T_solium", "H_microstoma"}


def distances(root):
    paths = {}
    def visit(node, path):
        path = dict(path)
        path[node] = node.length
        if not node.children:
            paths[node.name] = path
        for child in node.children:
            visit(child, path)
    visit(root, {})
    return {(a, b): sum(v for n, v in pa.items() if n not in pb)
                      + sum(v for n, v in pb.items() if n not in pa)
            for a, pa in paths.items() for b, pb in paths.items() if a < b}


class OrthologyPlotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.tree = self.root/"tree.nwk"
        self.tree.write_text(TREE)
        self.results = self.root/"results"
        folder = self.results/"Comparative_Genomics_Statistics"
        folder.mkdir(parents=True)
        self.species = list(m.LABELS)
        values = {sp: [1000, 900, 90, 2, 5] for sp in self.species}
        values["S_cotti"] = [9234, 8947, 96.9, 16, 692]
        values["S_solidus"] = [9796, 9540, 97.4, 20, 628]
        with (folder/"Statistics_PerSpecies.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([""]+self.species)
            for index, metric in enumerate(["Number of genes", "Number of genes in orthogroups",
                    "Percentage of genes in orthogroups", "Number of species-specific orthogroups"]):
                writer.writerow([metric]+[values[sp][index] for sp in self.species])
        with (folder/"Duplications_per_Species_Tree_Node.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Species Tree Node", "Duplications (50% support)"])
            for sp in self.species:
                writer.writerow([sp, values[sp][4]])

    def command(self, out):
        return [sys.executable, str(SCRIPT), "--tree", str(self.tree),
                "--results", str(self.results), "--out", str(out), "--font-scale", "1.5"]

    def test_binomial_labels_are_abbreviated_and_italic(self):
        for sp in self.species:
            if sp not in {"Spirometra_Aus1", "S_proliferum"}:
                parts = m.label_parts(sp)
                self.assertRegex(parts[0][0], r"^[A-Z][.] ")
                self.assertEqual(parts[0][1], "italic")
        self.assertEqual(m.label_parts("S_proliferum"), [("Sparganum proliferum", "italic")])

    def test_spirometra_rank_and_isolate_are_roman(self):
        self.assertEqual(m.label_parts("Spirometra_Aus1"),
                         [("Spirometra", "italic"), (" sp. Aus1", "normal")])

    def test_reroot_preserves_all_pairwise_distances(self):
        original = m.parse_newick(TREE)
        rooted, _, species = m.reroot_on_order_edge(original, CYCLO)
        self.assertEqual(species, set(self.species))
        before, after = distances(original), distances(rooted)
        for pair in before:
            self.assertAlmostEqual(before[pair], after[pair], places=12)

    def test_16_20_metrics_unchanged(self):
        metrics, dups = m.read_orthofinder_statistics(self.results, set(self.species))
        self.assertEqual((metrics["Number of species-specific orthogroups"]["S_cotti"],
                          metrics["Number of species-specific orthogroups"]["S_solidus"]), (16, 20))
        self.assertEqual((dups["S_cotti"], dups["S_solidus"]), (692, 628))

    def test_species_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, "species sets differ"):
            m.read_orthofinder_statistics(self.results, {"other"})

    def test_invalid_font_does_not_create_output(self):
        out = self.root/"invalid"
        command = self.command(out)+["--font-scale", "nan"]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(out.exists())

    def test_render_fonts_data_provenance_and_no_overwrite(self):
        out = self.root/"figure"
        sources = [self.tree]+list((self.results/"Comparative_Genomics_Statistics").iterdir())
        hashes = {p: m.sha256(p) for p in sources}
        command = self.command(out)
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertTrue((out/"orthology_summary_clean.pdf").read_bytes().startswith(b"%PDF"))
        audit = json.loads((out/"orthology_summary_clean.audit.json").read_text())
        self.assertEqual(audit["font_scale"], 1.5)
        self.assertEqual(audit["label_parts"]["Spirometra_Aus1"],
                         [["Spirometra", "italic"], [" sp. Aus1", "normal"]])
        texts = list(ET.parse(out/"orthology_summary_clean.svg").getroot().iter("{http://www.w3.org/2000/svg}text"))
        genus = next(t for t in texts if t.text == "Spirometra")
        suffix = next(t for t in texts if t.text == " sp. Aus1")
        self.assertIn("italic", genus.attrib["style"])
        self.assertNotIn("italic", suffix.attrib["style"])
        with (out/"orthology_summary_clean.source_data.tsv").open() as handle:
            rows = {r["species"]: r for r in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual((rows["S_cotti"]["assigned_proteins"], rows["S_solidus"]["assigned_proteins"]),
                         ("8947", "9540"))
        self.assertEqual(rows["S_proliferum"]["display_name"], "Sparganum proliferum")
        self.assertEqual(len(audit["branch_breaks"]), 3)
        self.assertTrue(audit["displayed_branches_compressed"])
        self.assertFalse(audit["source_tree_branch_lengths_modified"])
        self.assertEqual(audit["scale_applies_to"], "unbroken_segments_only")
        self.assertIn("not evolutionary", audit["tip_extensions"])
        self.assertEqual(hashes, {p: m.sha256(p) for p in sources})
        saved = m.sha256(out/"orthology_summary_clean.pdf")
        again = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(again.returncode, 0)
        self.assertEqual(saved, m.sha256(out/"orthology_summary_clean.pdf"))

    def test_optional_full_length_branches(self):
        out = self.root/"full_length"
        result = subprocess.run(self.command(out)+["--no-branch-breaks"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        audit = json.loads((out/"orthology_summary_clean.audit.json").read_text())
        self.assertEqual(audit["branch_breaks"], [])
        self.assertEqual(audit["scale_applies_to"], "all_branches")


if __name__ == "__main__":
    unittest.main()
