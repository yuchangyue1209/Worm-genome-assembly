#!/usr/bin/env python3
"""Prepare audited CAFE5 family-count files from OrthoFinder counts.

The analysis file retains families present in at least two selected species.
Families with at least 100 copies in any selected species are separated because
CAFE5 recommends estimating parameters on the smaller-family file. Explicitly
reviewed TE/artefact family IDs can be removed with --exclude-family-ids.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


def read_exclusions(paths: list[Path]) -> set[str]:
    excluded: set[str] = set()
    for path in paths:
        with path.open() as handle:
            for raw in handle:
                line = raw.strip()
                if line and not line.startswith("#"):
                    excluded.add(line.split()[0])
    return excluded


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-counts", type=Path, required=True,
                        help="OrthoFinder Orthogroups.GeneCount.tsv")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True,
                        help="Species to retain, in desired output order")
    parser.add_argument("--exclude-family-ids", type=Path, action="append", default=[])
    parser.add_argument("--large-family-threshold", type=int, default=100,
                        help="Separate a family if any selected species has this many copies")
    args = parser.parse_args()

    if args.large_family_threshold < 2:
        raise ValueError("--large-family-threshold must be at least 2")
    if len(args.species) != len(set(args.species)):
        raise ValueError("Duplicate species supplied")
    if args.out.exists():
        raise FileExistsError(f"Output already exists: {args.out}")

    excluded = read_exclusions(args.exclude_family_ids)
    args.out.mkdir(parents=True)

    small_path = args.out / "families_lt100.tsv"
    large_path = args.out / "families_ge100.tsv"
    audit_path = args.out / "family_filter_audit.tsv"
    summary = Counter()

    with args.gene_counts.open(newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup":
            raise ValueError("Expected an OrthoFinder Orthogroups.GeneCount.tsv file")
        missing = [sp for sp in args.species if sp not in reader.fieldnames]
        if missing:
            raise ValueError(f"Species absent from count table: {missing}")

        with small_path.open("x", newline="") as small_h, \
             large_path.open("x", newline="") as large_h, \
             audit_path.open("x", newline="") as audit_h:
            small = csv.writer(small_h, delimiter="\t", lineterminator="\n")
            large = csv.writer(large_h, delimiter="\t", lineterminator="\n")
            audit = csv.writer(audit_h, delimiter="\t", lineterminator="\n")
            header = ["Desc", "Family ID", *args.species]
            small.writerow(header)
            large.writerow(header)
            audit.writerow([
                "Family_ID", "Occupied_species", "Maximum_copy_number",
                "Selected_total", "Classification", "Reason",
            ])

            for row in reader:
                family = row["Orthogroup"]
                try:
                    counts = [int(row[sp]) for sp in args.species]
                except (TypeError, ValueError) as error:
                    raise ValueError(f"Non-integer count in {family}") from error
                if any(value < 0 for value in counts):
                    raise ValueError(f"Negative count in {family}")

                occupied = sum(value > 0 for value in counts)
                maximum = max(counts, default=0)
                total = sum(counts)
                record = ["(null)", family, *counts]

                if family in excluded:
                    classification, reason = "excluded", "reviewed_TE_or_artefact_family"
                elif occupied < 2:
                    classification, reason = "excluded", "fewer_than_two_selected_species"
                elif maximum >= args.large_family_threshold:
                    classification, reason = "large", f"copy_number_ge_{args.large_family_threshold}"
                    large.writerow(record)
                else:
                    classification, reason = "analysis", "retained_for_parameter_estimation"
                    small.writerow(record)

                summary[classification] += 1
                audit.writerow([family, occupied, maximum, total, classification, reason])

    unknown_exclusions = excluded - {
        row.split("\t", 1)[0]
        for row in audit_path.read_text().splitlines()[1:]
        if row
    }
    if unknown_exclusions:
        raise ValueError(
            "Excluded family IDs absent from input: " + ",".join(sorted(unknown_exclusions))
        )

    (args.out / "selected_species.txt").write_text("\n".join(args.species) + "\n")
    with (args.out / "summary.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Selected_species", len(args.species)])
        writer.writerow(["Families_for_parameter_estimation", summary["analysis"]])
        writer.writerow(["Large_families_separated", summary["large"]])
        writer.writerow(["Families_excluded", summary["excluded"]])
        writer.writerow(["Large_family_threshold", args.large_family_threshold])
        writer.writerow(["Explicit_excluded_family_IDs", len(excluded)])

    with (args.out / "sha256.tsv").open("x") as handle:
        for path in (small_path, large_path, audit_path):
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            handle.write(f"{digest}\t{path.name}\n")

    print(f"PASS: {small_path}")
    print(f"Species: {len(args.species)}")
    print(f"Families retained: {summary['analysis']}")
    print(f"Large families separated: {summary['large']}")
    print(f"Families excluded: {summary['excluded']}")


if __name__ == "__main__":
    main()
