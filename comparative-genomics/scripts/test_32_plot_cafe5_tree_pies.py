#!/usr/bin/env python3
"""Synthetic CAFE tree/count fixtures only; never substitutes for user data."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).with_name("32_plot_cafe5_tree_pies.py")
spec = importlib.util.spec_from_file_location("cafe_tree_plot", SCRIPT)
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)

DATED = "(((S_cotti:12,S_solidus:12):23,Spirometra_Aus1:35):25,(((E_granulosus:7,E_multilocularis:7):18,(T_multiceps:16.5,T_solium:16.5):8.5):15,H_microstoma:40):20);"
ASR = "(((S_cotti<1>_20:12,S_solidus<2>_20:12)<10>*_20:23,Spirometra_Aus1<3>_20:35)<11>_20:25,(((E_granulosus<4>_20:7,E_multilocularis<5>_20:7)<12>_20:18,(T_multiceps<6>_20:16.5,T_solium<7>_20:16.5)<13>_20:8.5)<14>_20:15,H_microstoma<8>_20:40)<15>_20:20)<9>_20;"


def fixture(directory, size=6):
    directory.mkdir(parents=True, exist_ok=True)
    dated, labelled = directory/"synthetic_dated.nwk", directory/"Gamma_asr.tre"
    dated.write_text(DATED)
    labelled.write_text("#NEXUS\nBEGIN TREES;\nTREE synthetic_family = [&R] "+ASR+"\nEND;\n")
    root, asr = m.parse_newick(DATED), m.first_asr_tree(labelled)
    m.map_nodes(root, asr)
    # Deliberately scramble columns and include root. Traversal order is not an ID map.
    nodes = sorted(m.walk(root), key=lambda n: (int(n.node_id)*7) % 16)
    names = {n.node_id: f"{next(iter(n.tips))}<{n.node_id}>" if not n.children else f"<{n.node_id}>" for n in nodes}
    columns = ["#FamilyID"]+[names[n.node_id] for n in nodes]
    changes, counts, family_tests = [], [], []
    for i in range(size):
        delta, count = {}, {}
        for node in m.walk(root):
            nid = int(node.node_id)
            value = 0 if node.parent is None else (1 if (i+nid) % 3 == 0 else -1 if (i+nid) % 3 == 1 else 0)
            delta[node.node_id] = value
            count[node.node_id] = 20 if node.parent is None else count[node.parent.node_id]+value
        fid = f"OG{i:07d}"
        changes.append([fid]+[delta[n.node_id] for n in nodes])
        counts.append([fid]+[count[n.node_id] for n in nodes])
        family_tests.append([fid, ["0", "0.001", "0.01", "0.04", "0.5", "NA"][i % 6]])
    for stem, rows, header in [("Gamma_change.tab", changes, columns),
                               ("Gamma_count.tab", counts, columns),
                               ("Gamma_family_results.txt", family_tests, ["#FamilyID", "pvalue"])]:
        with (directory/stem).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(header)
            writer.writerows(rows)
    stats, _, _ = m.calculate_stats(root, directory/"Gamma_change.tab", directory/"Gamma_count.tab")
    with (directory/"Gamma_clade_results.txt").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["#Taxon_ID", "Increase", "Decrease"])
        for row in stats.values():
            writer.writerow([row["Branch"], row["Expanded"], row["Contracted"]])
    return dated, labelled


class CafePlotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dated, self.asr = fixture(self.root/"results")
        self.results = self.dated.parent

    def mapped(self):
        root = m.parse_newick(self.dated.read_text())
        height = m.map_nodes(root, m.first_asr_tree(self.asr))
        return root, height

    def stats(self, threshold=None):
        root, _ = self.mapped()
        return m.calculate_stats(root, self.results/"Gamma_change.tab", self.results/"Gamma_count.tab",
                                 self.results/"Gamma_family_results.txt", threshold)

    def test_explicit_ids_and_roundtrip_clades(self):
        root, height = self.mapped()
        self.assertEqual(height, 60)
        self.assertEqual(root.node_id, "9")
        pair = next(n for n in m.walk(root) if n.tips == frozenset(["S_cotti", "S_solidus"]))
        self.assertEqual(pair.node_id, "10")
        self.assertEqual(height-pair.depth, 12)
        leaves = m.layout(root, m.TIP_ORDER)
        self.assertEqual([next(iter(n.tips)) for n in leaves], m.TIP_ORDER)

    def test_color_vision_friendly_palette(self):
        self.assertEqual(m.COLORS, {"Expanded": "#0072B2", "Contracted": "#D55E00", "Unchanged": "#D9D9D9"})
        self.assertEqual(m.PALETTE_SOURCE, "https://jfly.uni-koeln.de/color/")

    def test_wrong_rooted_topology_rejected(self):
        altered = DATED.replace("S_solidus", "TEMP").replace("Spirometra_Aus1", "S_solidus").replace("TEMP", "Spirometra_Aus1")
        with self.assertRaisesRegex(ValueError, "topology"):
            m.map_nodes(m.parse_newick(altered), m.first_asr_tree(self.asr))

    def test_wrong_calibration_scale_rejected(self):
        altered = DATED.replace(":12,", ":13,")
        with self.assertRaisesRegex(ValueError, "branch length"):
            m.map_nodes(m.parse_newick(altered), m.first_asr_tree(self.asr))

    def test_family_counts_not_asr_one_family_copies(self):
        stats, total, selected = self.stats()
        self.assertEqual((total, selected), (6, 6))
        self.assertEqual(len(stats), 14)
        self.assertNotIn("9", stats)
        for row in stats.values():
            self.assertEqual((row["Expanded"], row["Contracted"], row["Unchanged"]), (2, 2, 2))

    def test_nominal_subset_strict_boundary_and_denominator(self):
        stats, total, selected = self.stats(.01)
        self.assertEqual((total, selected), (6, 2))
        for row in stats.values():
            self.assertEqual(row["Expanded"]+row["Contracted"]+row["Unchanged"], 2)
            self.assertEqual(row["All_expanded"], 2)

    def test_clade_summary_crosscheck(self):
        stats, _, _ = self.stats()
        path = self.results/"Gamma_clade_results.txt"
        m.check_clade_summary(path, stats, "9")
        path.write_text(path.read_text().replace("\t2\t2", "\t3\t2", 1))
        with self.assertRaisesRegex(ValueError, "disagrees"):
            m.check_clade_summary(path, stats, "9")

    def test_wrong_tip_node_label_rejected(self):
        path = self.results/"Gamma_change.tab"
        path.write_text(path.read_text().replace("S_cotti<1>", "S_solidus<1>"))
        with self.assertRaisesRegex(ValueError, "species/node"):
            self.stats()

    def test_parent_count_consistency_is_checked(self):
        path = self.results/"Gamma_change.tab"
        header, rows = m.read_table(path)
        next(iter(rows.values()))["S_cotti<1>"] = "99"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, header, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows.values())
        with self.assertRaisesRegex(ValueError, "parent reconstruction mismatch"):
            self.stats()

    def test_unknown_clade_ids_are_not_silently_ignored(self):
        stats, _, _ = self.stats()
        path = self.results/"Gamma_clade_results.txt"
        with path.open("a") as handle:
            handle.write("<999>\t0\t0\n")
        with self.assertRaisesRegex(ValueError, "Unknown clade-summary"):
            m.check_clade_summary(path, stats, "9")

    def test_missing_asr_does_not_guess_mapping_or_create_output(self):
        self.asr.unlink()
        out = self.root/"figures"
        result = subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
            "--tree", str(self.dated), "--out", str(out)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("node-labelled ASR tree required", result.stderr)
        self.assertFalse(out.exists())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib unavailable")
    def test_full_plot_and_source_preservation(self):
        out = self.root/"figures"
        sources = {p: m.sha256(p) for p in self.results.iterdir()}
        command = [sys.executable, str(SCRIPT), "--results", str(self.results),
            "--tree", str(self.dated), "--out", str(out), "--dpi", "100", "--font-scale", "1.15"]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertTrue((out/"cafe5_tree_pies.pdf").read_bytes().startswith(b"%PDF"))
        self.assertTrue((out/"cafe5_tree_pies.png").read_bytes().startswith(b"\x89PNG"))
        self.assertIn("Million years ago", (out/"cafe5_tree_pies.svg").read_text())
        _, rows = m.read_table(out/"branch_family_counts.tsv")
        self.assertEqual(len(rows), 14)
        config = json.loads((out/"figure_provenance.json").read_text())
        self.assertEqual(config["all_reconstructed_families"], 6)
        self.assertFalse(config["branch_probabilities_used_as_pvalues"])
        self.assertEqual(config["clade_summary_crosscheck"], "passed")
        self.assertEqual(config["palette"], m.COLORS)
        self.assertIn("signed +/-", config["redundant_encoding"])
        svg = (out/"cafe5_tree_pies.svg").read_text().lower()
        for color in m.COLORS.values():
            self.assertIn(color.lower(), svg)
        self.assertEqual(len(config["internal_annotation_offsets"]), 6)
        self.assertEqual(sources, {p: m.sha256(p) for p in self.results.iterdir()})
        again = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(again.returncode, 0)
        self.assertIn("Output exists", again.stderr)


if __name__ == "__main__":
    unittest.main()
