#!/usr/bin/env python3
"""Synthetic tests of family-wide screening; no actual CAFE jobs are launched."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).with_name("29_summarize_cafe5.py")
spec = importlib.util.spec_from_file_location("summary", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class SummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.results = self.root/"results"
        self.results.mkdir()
        self.input = self.root/"input.tsv"
        self.families = [
            ("zero", "0", 1, -1),
            ("small", "0.001", -1, 1),
            ("boundary", "0.01", 1, 0),
            ("marginal", "0.04", 1, 1),
            ("unchanged", "0.0001", 0, 0),
            ("not_tested", "NA", 1, -1),
            ("high", "0.8", 1, 0),
        ]
        (self.results/"Gamma_family_results.txt").write_text(
            "#FamilyID\tpvalue\tSignificant at 0.05\n"+"".join(
                f"{og}\t{p}\t{'y' if p != 'NA' and float(p) < .05 else 'n'}\n"
                for og,p,_,_ in self.families))
        header = "#FamilyID\tS_cotti<1>\tS_solidus<2>\t<3>\n"
        (self.results/"Gamma_change.tab").write_text(header+"".join(
            f"{og}\t{c}\t{s}\t0\n" for og,_,c,s in self.families))
        (self.results/"Gamma_count.tab").write_text(header+"".join(
            f"{og}\t{1+c}\t{1+s}\t1\n" for og,_,c,s in self.families))
        # A large posterior, a tiny posterior, and absent values must NOT gate candidates.
        (self.results/"Gamma_branch_probabilities.tab").write_text(
            header+"zero\t0.99\t0.8\tNA\nboundary\t0.00001\tNA\tNA\n")
        self.input.write_text("Desc\tFamily ID\tS_cotti\tS_solidus\n"+"".join(
            f"(null)\t{og}\t{1+c}\t{1+s}\n" for og,_,c,s in self.families)+
            "(null)\troot_filtered\t1\t0\n")
        self.n = 0

    def read(self, path):
        with path.open(newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def run_summary(self, *options, success=True, out=None):
        self.n += 1
        out = out or self.root/f"summary{self.n}"
        cmd = [sys.executable, str(SCRIPT), "--results", str(self.results),
               "--families", str(self.input), "--out", str(out), *options]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return out, result

    def test_nominal_default_ignores_branch_probabilities(self):
        out, _ = self.run_summary()
        rows = self.read(out/"focal_nominal_candidates.tsv")
        for sp in ["S_cotti", "S_solidus"]:
            self.assertEqual({r["Orthogroup"] for r in rows if r["Species"] == sp}, {"zero", "small"})
        zero = next(r for r in rows if r["Species"] == "S_cotti" and r["Orthogroup"] == "zero")
        self.assertEqual(zero["Branch_probability_reported"], "0.99")
        self.assertNotIn("Branch_pvalue_raw", zero)
        self.assertNotIn("Branch_p_pass", zero)
        self.assertEqual(zero["Family_test_scope"], "family_wide_not_branch_specific")
        summary = self.read(out/"focal_summary.tsv")
        self.assertEqual(summary[0]["Nominal_family_p_expanded"], "1")
        self.assertEqual(summary[0]["Nominal_family_p_contracted"], "1")
        self.assertEqual(summary[0]["Tested_families"], "6")
        self.assertEqual(summary[0]["BH_floor_sensitivity_expanded"], "NA")

    def test_optional_branch_table(self):
        (self.results/"Gamma_branch_probabilities.tab").unlink()
        out, _ = self.run_summary()
        rows = self.read(out/"focal_nominal_candidates.tsv")
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(r["Branch_probability_reported"] == "NA" for r in rows))

    def test_threshold_strict_boundary_and_reanalysis(self):
        out, _ = self.run_summary("--family-p-threshold", "0.05")
        rows = self.read(out/"focal_nominal_candidates.tsv")
        cotti = {r["Orthogroup"] for r in rows if r["Species"] == "S_cotti"}
        self.assertEqual(cotti, {"zero", "small", "boundary", "marginal"})
        self.assertFalse({"unchanged", "not_tested", "high"} & cotti)

    def test_fdr_threshold_separate_from_nominal_threshold(self):
        out, _ = self.run_summary("--fdr-threshold", "0.001")
        candidates = self.read(out/"focal_nominal_candidates.tsv")
        self.assertEqual(len(candidates), 4)
        small = next(r for r in candidates if r["Orthogroup"] == "small")
        self.assertEqual(small["Family_p_pass"], "1")
        self.assertEqual(small["Family_BH_reported_pass"], "0")
        self.assertAlmostEqual(float(small["Family_qvalue_BH_reported"]), .002)

    def test_floor_sensitivity_does_not_replace_raw_p_or_candidates(self):
        out, _ = self.run_summary("--pvalue-floor", "0.02")
        candidates = self.read(out/"focal_nominal_candidates.tsv")
        self.assertEqual(len(candidates), 4)
        row = next(r for r in candidates if r["Orthogroup"] == "zero")
        self.assertEqual(row["Family_pvalue"], "0.0")
        self.assertEqual(row["Family_qvalue_BH_reported"], "0.0")
        self.assertEqual(row["Family_pvalue_floor_sensitivity"], "0.02")
        self.assertGreater(float(row["Family_qvalue_BH_floor_sensitivity"]), 0)
        settings = json.loads((out/"summary_settings.json").read_text())
        self.assertFalse(settings["branch_probabilities_used_as_pvalues"])
        self.assertEqual(settings["family_p_threshold"], .01)
        self.assertEqual(settings["bh_q_threshold"], .05)

    def test_omissions_and_source_protection(self):
        inputs = [*self.results.iterdir(), self.input]
        hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
        out, _ = self.run_summary()
        omitted = self.read(out/"input_families_not_in_output.tsv")
        self.assertEqual([r["Orthogroup"] for r in omitted], ["root_filtered"])
        self.assertEqual(hashes, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
        self.run_summary(out=out, success=False)
        self.assertTrue((out/"focal_summary.tsv").is_file())

    def test_bh_known_values(self):
        q = m.bh_adjust(dict(a=.01, b=.04, c=.03, d=.002))
        for key, value in dict(a=.02, b=.04, c=.04, d=.008).items():
            self.assertAlmostEqual(q[key], value)

    def test_invalid_thresholds(self):
        for option, value in [("--family-p-threshold", "0"), ("--family-p-threshold", "1"),
                              ("--fdr-threshold", "nan"), ("--pvalue-floor", "nan")]:
            self.run_summary(option, value, success=False)

    def test_unknown_branch_family_is_rejected(self):
        with (self.results/"Gamma_branch_probabilities.tab").open("a") as f:
            f.write("unknown\t0.5\t0.5\tNA\n")
        _, result = self.run_summary(success=False)
        self.assertIn("unknown family IDs", result.stderr)


if __name__ == "__main__":
    unittest.main()
