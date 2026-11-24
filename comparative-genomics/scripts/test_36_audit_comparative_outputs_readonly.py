#!/usr/bin/env python3
"""Small read-only inventory tests with synthetic files, never server data."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("36_audit_comparative_outputs_readonly.py")


class InventoryTests(unittest.TestCase):
    def test_inventory_preserves_files_and_distinguishes_presence(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            code = root/"comparative-genomics/scripts/example.py"
            code.parent.mkdir(parents=True)
            code.write_text("# test fixture only\n")
            manifest = root/"manifest.json"
            manifest.write_text(json.dumps({"files": [
                {"path": str(code.relative_to(root)), "sha256": hashlib.sha256(code.read_bytes()).hexdigest()},
                {"path": "comparative-genomics/scripts/not_present.py", "sha256": "0"*64}]}))
            results = root/"07_gene_family/cafe5/runs/screening8_16.5Ma_k3_v2/02_gamma_k3"
            results.mkdir(parents=True)
            (results/"Gamma_results.txt").write_text(
                "Model Gamma Final Likelihood (-lnL): 123.4\nLambda: 0.002\nAlpha: 0.5\n")
            (results/"Gamma_change.tab").write_text("# synthetic presence fixture\n")
            paths = [p for p in root.rglob("*") if p.is_file()]
            before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
            proc = subprocess.run([sys.executable, str(SCRIPT), "--root", str(root),
                                   "--code-manifest", str(manifest)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("PRESENT_NOT_VALIDATED", proc.stdout)
            self.assertIn("PARSED_NOT_CONVERGENCE_PROOF", proc.stdout)
            self.assertIn("MISSING_OR_EMPTY", proc.stdout)
            self.assertIn("not_present.py\tMISSING", proc.stdout)
            self.assertIn("READ_ONLY", proc.stdout)
            after = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in root.rglob("*") if p.is_file()}
            self.assertEqual(before, after)

    def test_manifest_cannot_escape_project(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            manifest = root/"manifest.json"
            manifest.write_text(json.dumps({"files": [{"path": "../outside.py", "sha256": "0"*64}]}))
            proc = subprocess.run([sys.executable, str(SCRIPT), "--root", str(root),
                                   "--code-manifest", str(manifest)], capture_output=True, text=True)
            self.assertIn("CHECK_FAILED", proc.stdout)
            self.assertIn("relative to the project", proc.stdout)


if __name__ == "__main__":
    unittest.main()
