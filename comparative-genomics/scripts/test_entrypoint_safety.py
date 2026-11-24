#!/usr/bin/env python3
"""Synthetic safety tests: never invoke real BUSCO or OrthoFinder."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parent

class EntryPointSafetyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="comparative-entry-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.primary = self.root / "std" / "proteins_primary"
        self.primary.mkdir(parents=True)
        for species in ("A", "B"):
            (self.primary / (species + ".faa")).write_text(">p1\nMAAA\n")
        self.config = self.root / "config.sh"
        self.write_config("A B")
        self.mock("orthofinder", '#!/bin/bash\nprintf "%s\\n" "$@" > "$FIXTURE_ROOT/orthofinder.args"\n')
        self.mock("sha256sum", "#!" + sys.executable + "\nimport hashlib,sys\nfrom pathlib import Path\nfor arg in sys.argv[1:]:\n print(hashlib.sha256(Path(arg).read_bytes()).hexdigest()+'  '+arg)\n")
        self.env = os.environ.copy()
        for key in ("INCLUDE_SCOTTI_RELEASE", "INCLUDE_SCOTTI_CANDIDATE", "BUSCO_LINEAGE_OVERRIDE"):
            self.env.pop(key, None)
        self.env.update(CONFIG=str(self.config), FIXTURE_ROOT=str(self.root),
                        ORTHO_RUN_ID="fixture", PATH=str(self.bin) + os.pathsep + str(Path(sys.executable).parent) + os.pathsep + self.env["PATH"])

    def mock(self, name, content):
        target = self.bin / name
        target.write_text(content)
        target.chmod(0o755)

    def write_config(self, species):
        self.config.write_text(
            'PROJECT="' + str(self.root) + '"\nSTD="$PROJECT/std"\n'
            'ORTHO_ROOT="$PROJECT/orthofinder"\nQC="$PROJECT/qc"\n'
            'THREADS=1\nBUSCO_LINEAGE=metazoa_odb12\nSPECIES=(' + species + ')\n')

    def run_script(self, name):
        return subprocess.run(["bash", str(SCRIPTS / name)], env=self.env,
                              text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def test_run_freezes_inputs_and_refuses_reuse(self):
        result = self.run_script("03_orthofinder.sh")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        frozen = self.root / "orthofinder/runs/fixture/input/A.faa"
        self.assertFalse(frozen.is_symlink())
        before = frozen.read_bytes()
        (self.primary / "A.faa").write_text(">changed\nMCCC\n")
        self.assertEqual(frozen.read_bytes(), before)
        args = (self.root / "orthofinder.args").read_bytes()
        result = self.run_script("03_orthofinder.sh")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Run exists", result.stdout)
        self.assertEqual((self.root / "orthofinder.args").read_bytes(), args)
        self.assertEqual(frozen.read_bytes(), before)

    def test_duplicate_ids_stop_before_job_or_run_directory(self):
        (self.primary / "A.faa").write_text(">same\nMA\n>same\nMC\n")
        result = self.run_script("03_orthofinder.sh")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate protein ID", result.stderr)
        self.assertFalse((self.root / "orthofinder.args").exists())
        self.assertFalse((self.root / "orthofinder/runs/fixture").exists())

    def test_wrong_cotti_baseline_stops_before_job(self):
        self.write_config("A S_cotti")
        (self.primary / "S_cotti.faa").write_text(
            "".join(">g" + str(i) + "\nMA\n" for i in range(9234)))
        result = self.run_script("03_orthofinder.sh")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("differs from the selected 9234-protein release", result.stderr)
        self.assertFalse((self.root / "orthofinder.args").exists())
        self.assertFalse((self.root / "orthofinder/runs/fixture").exists())

    def test_busco_lineages_have_separate_output_directories(self):
        self.mock("busco", '#!/bin/bash\nprintf "%s\\n" "$@" >> "$FIXTURE_ROOT/busco.args"\n')
        for lineage in ("metazoa_odb10", "metazoa_odb12"):
            self.env["BUSCO_LINEAGE_OVERRIDE"] = lineage
            result = self.run_script("02b_busco_proteins.sh")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue((self.root / "qc/busco_protein" / lineage).is_dir())
        args = (self.root / "busco.args").read_text()
        self.assertIn(str(self.root / "qc/busco_protein/metazoa_odb10"), args)
        self.assertIn(str(self.root / "qc/busco_protein/metazoa_odb12"), args)

if __name__ == "__main__":
    unittest.main()
