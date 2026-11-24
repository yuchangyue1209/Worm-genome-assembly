#!/usr/bin/env python3
"""Inventory known comparative outputs and code versions without running analyses.

Standard library only. All reports go to stdout. No files, process state,
CAFE STATUS flags, FASTAs, links or result directories are modified.
Presence and parsable summaries are not scientific completion validation.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re


ROOT = Path("/work/cyu/Schistocephalus_comparative_genomics")
BASE = "04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/OrthoFinder/Results_Nov09"
SENS = "04_orthofinder/runs/exploratory11_scotti_without32_v2/input/OrthoFinder/Results_Nov12"
CAFE = "07_gene_family/cafe5"
FIT = CAFE+"/runs/screening8_16.5Ma_k3_v2/02_gamma_k3"
SUMMARY = CAFE+"/summaries/screening8_16.5Ma_k3_v2_focal_nominal_p0.01_schema2_v2"
TABLE = "03_qc/release_table1_11species_v5_exonfix"
COTTI_SHA = "07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371"
PROBES = [
    ("BUSCO", "22-assessment combined table", ("03_qc/final_genome_busco_odb12_v1/eleven_species.genome_and_protein.odb12.tsv",)),
    ("BUSCO", "batch completion record", ("03_qc/final_genome_busco_odb12_v1/COMPLETE.json",)),
    ("Table1", "preserved v4 table", ("03_qc/release_table1_11species_v4/Table1_11species.tsv",)),
    ("Table1", "corrected v5 table", (TABLE+"/Table1_11species.tsv",)),
    ("Table1", "corrected exon audit", (TABLE+"/exon_transcript_audit.tsv",)),
    ("Table1", "v5 collection settings", (TABLE+"/collection.json",)),
    ("OrthoFinder", "baseline membership", (BASE+"/Orthogroups/Orthogroups.tsv",)),
    ("OrthoFinder", "without32 v2 membership", (SENS+"/Orthogroups/Orthogroups.tsv",)),
    ("Phylogeny", "without32 v2 IQ-TREE", ("04_orthofinder/phylogeny_without32_v2/concatenated.treefile",)),
    ("Functional", "eleven-species coverage table", ("05_functional_annotation/annotation_completeness/functional_coverage.tsv",)),
    ("Functional", "proliferum InterPro", ("05_functional_annotation/interpro/S_proliferum.tsv",)),
    ("Functional", "proliferum eggNOG", ("05_functional_annotation/eggnog/S_proliferum.emapper.annotations",)),
    ("Species-specific", "baseline focal OG membership", ("07_gene_family/species_specific_baseline9234/focal_species_specific_orthogroups.tsv",)),
    ("Species-specific", "15-protein RNA support", ("07_gene_family/species_specific_baseline9234/rna_repeat_validation/priority_candidate_rnaseq_support.tsv",)),
    ("CAFE", "eight-species family matrix", (CAFE+"/screening8_input_v1/families_lt100.tsv",)),
    ("CAFE", "dated screening tree", (CAFE+"/time_trees/screening8_dated_16.5Ma.nwk",)),
    ("CAFE", "Gamma family tests", (FIT+"/Gamma_family_results.txt",)),
    ("CAFE", "Gamma changes", (FIT+"/Gamma_change.tab", FIT+"/Gamma_change.txt")),
    ("CAFE", "Gamma reconstructed counts", (FIT+"/Gamma_count.tab", FIT+"/Gamma_count.txt")),
    ("CAFE", "Gamma ASR node mapping", (FIT+"/Gamma_asr.tre",)),
    ("CAFE", "new schema2 summary", (SUMMARY+"/focal_summary.tsv",)),
    ("CAFE", "new schema2 settings", (SUMMARY+"/summary_settings.json",)),
    ("CAFE-review", "post-fix priority4 summary", (CAFE+"/candidate_review/priority4_baseline9234_v2/family_review_summary.tsv",)),
    ("CAFE-review", "post-fix sequence audit", (CAFE+"/candidate_review/priority4_baseline9234_v2/protein_sequence_match_audit.tsv",)),
    ("Figure2", "all-change PDF", (CAFE+"/figures/screening8_16.5Ma_k3_v2_all_changes_cvd_v1/cafe5_tree_pies.pdf",)),
    ("Figure2", "figure provenance", (CAFE+"/figures/screening8_16.5Ma_k3_v2_all_changes_cvd_v1/figure_provenance.json",)),
    ("Secretome", "integrated summary", ("06_secretome/integrated_final/secretome_summary.tsv",)),
    ("MEROPS", "classical-secretome summary", ("06_secretome/merops/merops_confirmed_secretome_summary_all_species.tsv",)),
]


def report(section, item, status, detail):
    print("\t".join(str(v).replace("\t", " ").replace("\n", " ")
                    for v in (section, item, status, detail)))


def code_audit(root, manifest):
    if not manifest.is_file():
        report("Code", "hash manifest", "NOT_AVAILABLE", manifest)
        return
    try:
        payload = json.loads(manifest.read_text())
        counts = Counter()
        for row in payload["files"]:
            relative = Path(row["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Manifest paths must be relative to the project")
            path = root/relative
            if not path.is_file():
                state = "MISSING"
            elif path.stat().st_size > 12*1024*1024:
                state = "TOO_LARGE_FOR_CODE_AUDIT"
            else:
                state = ("MATCH" if hashlib.sha256(path.read_bytes()).hexdigest()
                         == row["sha256"] else "DIFFERENT")
            counts[state] += 1
            if state != "MATCH":
                report("Code", relative, state, "Compared with local synchronization manifest")
        report("Code", "hash comparison", "INVENTORY_ONLY", dict(counts))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        report("Code", "hash manifest", "CHECK_FAILED", exc)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--code-manifest", type=Path)
    args = parser.parse_args()
    root = args.root
    print("Section\tItem\tStatus\tDetails")
    code_audit(root, args.code_manifest or root/"comparative-genomics/data/analysis_sync_manifest.json")
    for section, name, options in PROBES:
        paths = [root/p for p in options]
        found = next((p for p in paths if p.is_file() and p.stat().st_size), None)
        report(section, name, "PRESENT_NOT_VALIDATED" if found else "MISSING_OR_EMPTY",
               f"{found} ({found.stat().st_size} bytes)" if found else " | ".join(map(str, paths)))
    for label, relative in (("baseline", BASE), ("without32_v2", SENS)):
        path = root/relative/"Orthogroups/Orthogroups.tsv"
        if path.is_file():
            try:
                with path.open() as handle:
                    fields = next(csv.reader(handle, delimiter="\t"))
                species = fields[1:]
                good = len(species) == 11 and len(set(species)) == 11 and "S_cotti" in species
                report("OrthoFinder", label+" header", "11_SPECIES_HEADER" if good else "CHECK_FAILED",
                       ",".join(species))
            except (OSError, StopIteration) as exc:
                report("OrthoFinder", label+" header", "CHECK_FAILED", exc)
    for version in (2, 3, 4):
        path = root/CAFE/f"runs/screening8_16.5Ma_k3_v{version}/02_gamma_k3/Gamma_results.txt"
        if not path.is_file():
            report("CAFE-fit", f"v{version}", "MISSING", path)
            continue
        try:
            text = path.read_text()
            values = {}
            for label, pattern in (
                    ("nll", r"^Model Gamma Final Likelihood \(-lnL\):\s*(\S+)"),
                    ("lambda", r"^Lambda:\s*(\S+)"), ("alpha", r"^Alpha:\s*(\S+)")):
                match = re.search(pattern, text, re.MULTILINE)
                if not match:
                    raise ValueError("Missing "+label)
                values[label] = float(match.group(1))
            if not all(math.isfinite(v) for v in values.values()) or min(values["lambda"], values["alpha"]) <= 0:
                raise ValueError("Non-finite or invalid fit parameters")
            report("CAFE-fit", f"v{version}", "PARSED_NOT_CONVERGENCE_PROOF", values)
        except (OSError, ValueError) as exc:
            report("CAFE-fit", f"v{version}", "CHECK_FAILED", exc)
    for label, relative in (("reported_nominal_v1", CAFE+"/summaries/screening8_16.5Ma_k3_v2_focal_nominal_p0.01_v1"),
                             ("new_schema2_v2", SUMMARY)):
        path = root/relative/"summary_settings.json"
        if path.is_file():
            try:
                settings = json.loads(path.read_text())
                good = settings.get("schema_version") == 2 and settings.get("branch_probabilities_used_as_pvalues") is False
                report("CAFE-summary", label, "SCHEMA2" if good else "OTHER_OR_LEGACY_SCHEMA", settings)
            except (OSError, ValueError) as exc:
                report("CAFE-summary", label, "CHECK_FAILED", exc)
    for label, path in (
            ("frozen_OrthoFinder_cotti", root/"04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/S_cotti.faa"),
            ("shared_release_cotti", Path("/work/shared/Scotti_release_v1.0_preNCBI_20261118/proteins/Scotti_proteins.longest_isoform.v1.0.fa"))):
        if path.is_file():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            report("Baseline", label, "MATCH_RECORDED_SHA256" if digest == COTTI_SHA else "DIFFERENT_FROM_RECORDED_BASELINE", digest)
    report("Audit", "scope", "READ_ONLY", "No jobs launched; no STATUS flags or files modified. Presence alone does not prove completion or eleven-species coverage.")


if __name__ == "__main__":
    main()
