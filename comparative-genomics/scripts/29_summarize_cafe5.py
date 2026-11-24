#!/usr/bin/env python3
"""Screen family-wide CAFE5 tests and describe reconstructed focal copy changes.

Family BH q-values cover every tested family in the supplied model output.
Branch probabilities are retained for reference, NEVER treated as p-values or
used as a candidate filter. In particular, Gamma branch outputs are posterior
probabilities, not branch-specific p-values (CAFE5 developer discussion #232).
An optional p-value floor provides a labelled resolution-sensitivity scenario,
not a replacement estimate of the true empirical p-values.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re


def read_table(path: Path):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        while header and not header[-1].strip():
            header.pop()
        header = [name.strip().lstrip("#") for name in header]
        if not header or len(header) != len(set(header)):
            raise ValueError(f"Missing or duplicate header fields: {path}")
        rows = {}
        for line_number, cells in enumerate(reader, 2):
            if not cells or not any(cell.strip() for cell in cells):
                continue
            while len(cells) > len(header) and not cells[-1].strip():
                cells.pop()
            if len(cells) != len(header):
                raise ValueError(f"Unexpected field count: {path}:{line_number}")
            cells = [cell.strip() for cell in cells]
            if not cells[0] or cells[0] in rows:
                raise ValueError(f"Empty/duplicate family identifier: {path}:{line_number}")
            rows[cells[0]] = dict(zip(header, cells))
    return header, rows


def resolve_table(directory: Path, stem: str) -> Path:
    for extension in ("tab", "txt", "tsv"):
        path = directory / f"{stem}.{extension}"
        if path.is_file() and path.stat().st_size:
            return path
    raise FileNotFoundError(f"Missing nonempty {directory}/{stem}.tab/.txt/.tsv")


def probability(raw: str) -> float | None:
    if raw.upper() in {"NA", "N/A", "", "NAN", "-NAN"}:
        return None
    value = float(raw)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"Invalid probability: {raw}")
    return value


def bh_adjust(values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(values, key=lambda key: (values[key], key))
    size = len(ordered)
    adjusted = {}
    running_minimum = 1.0
    for rank in range(size, 0, -1):
        family = ordered[rank - 1]
        running_minimum = min(running_minimum, values[family] * size / rank)
        adjusted[family] = running_minimum
    return adjusted


def write_table(path: Path, columns: list[str], rows: list[dict]):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else row[key]
                             for key in columns})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True,
                        help="Directory containing Gamma_results.txt and reconstruction tables")
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--species", nargs="+", default=["S_cotti", "S_solidus"])
    parser.add_argument("--model", choices=["Gamma", "Base"], default="Gamma")
    parser.add_argument("--families", type=Path,
                        help="Original CAFE input, for auditing families omitted internally")
    parser.add_argument("--pvalue-floor", type=float,
                        help="Sensitivity only: floor reported family p-values before BH; retain raw values")
    parser.add_argument("--family-p-threshold", "--threshold", dest="family_p_threshold",
                        type=float, default=0.01,
                        help="Nominal family-wide p threshold; default 0.01. Not a branch-specific test.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05,
                        help="Separate BH q threshold for supplementary diagnostics; default 0.05")
    args = parser.parse_args()
    if not 0 < args.family_p_threshold < 1:
        raise ValueError("--family-p-threshold must be between zero and one")
    if not 0 < args.fdr_threshold < 1:
        raise ValueError("--fdr-threshold must be between zero and one")
    if args.pvalue_floor is not None and not 0 < args.pvalue_floor < 1:
        raise ValueError("--pvalue-floor must be finite and between zero and one")
    if len(set(args.species)) != len(args.species):
        raise ValueError("Duplicate focal species")
    if args.out.exists():
        raise FileExistsError(args.out)

    paths = {kind: resolve_table(args.results, f"{args.model}_{kind}")
             for kind in ("family_results", "change", "count")}
    family_header, families = read_table(paths["family_results"])
    change_header, changes = read_table(paths["change"])
    count_header, counts = read_table(paths["count"])
    # The nominal screen does not require a branch-probability table.
    branch_header, branches = [], {}
    try:
        paths["branch_probabilities"] = resolve_table(args.results, f"{args.model}_branch_probabilities")
    except FileNotFoundError:
        pass
    else:
        branch_header, branches = read_table(paths["branch_probabilities"])
    p_column = next((name for name in family_header if name.lower() == "pvalue"), None)
    if p_column is None:
        raise ValueError(f"Cannot find pvalue column: {family_header}")
    if not families or set(families) != set(changes) or set(families) != set(counts):
        raise ValueError("Family results, count and change tables must contain identical families")
    if set(branches) - set(families):
        raise ValueError("Branch probabilities include unknown family IDs")

    pvalues = {family: probability(row[p_column]) for family, row in families.items()}
    tested = {family: value for family, value in pvalues.items() if value is not None}
    zero_pvalues = sum(value == 0 for value in tested.values())
    positive_pvalues = [value for value in tested.values() if value > 0]
    reported_qvalues = bh_adjust(tested)
    sensitivity_pvalues = ({family: max(value, args.pvalue_floor)
                           for family, value in tested.items()}
                          if args.pvalue_floor is not None else {})
    sensitivity_qvalues = bh_adjust(sensitivity_pvalues)
    ordered = sorted(families, key=lambda family: (
        pvalues[family] if pvalues[family] is not None else 2, family))
    focal_rows = []
    summary_rows = []
    for species in args.species:
        matches = [name for name in change_header[1:]
                   if re.fullmatch(re.escape(species) + r"(?:<\d+>)?", name)]
        if len(matches) != 1:
            raise ValueError(f"Expected one branch for {species}; found {matches}")
        column = matches[0]
        if column not in count_header or (branch_header and column not in branch_header):
            raise ValueError(f"Branch label absent from count/probability tables: {column}")
        records = []
        for family in ordered:
            delta = int(changes[family][column])
            observed = int(counts[family][column])
            parent = observed - delta
            if min(observed, parent) < 0:
                raise ValueError(f"Negative observed/parent count: {species}/{family}")
            pvalue = pvalues[family]
            reported_q = reported_qvalues.get(family)
            sensitivity_q = sensitivity_qvalues.get(family)
            branch_probability = probability(branches.get(family, {}).get(column, "NA"))
            family_raw = pvalue is not None and pvalue < args.family_p_threshold
            reported_bh_pass = reported_q is not None and reported_q < args.fdr_threshold
            sensitivity_bh_pass = (sensitivity_q < args.fdr_threshold
                                   if sensitivity_q is not None else None)
            records.append(dict(
                Species=species, Orthogroup=family, Branch=column,
                Observed_copies=observed, Parent_reconstructed_copies=parent,
                Copy_change=delta,
                Direction="Expansion" if delta > 0 else "Contraction" if delta < 0 else "Unchanged",
                Family_pvalue=pvalue,
                Family_qvalue_BH_reported=reported_q,
                Family_pvalue_floor_sensitivity=sensitivity_pvalues.get(family),
                Family_qvalue_BH_floor_sensitivity=sensitivity_q,
                Branch_probability_reported=branch_probability,
                Family_p_pass=int(family_raw),
                Family_BH_reported_pass=int(reported_bh_pass),
                Family_BH_floor_sensitivity_pass=(int(sensitivity_bh_pass)
                                                if sensitivity_bh_pass is not None else None),
                Nominal_family_change_candidate=int(delta != 0 and family_raw),
                BH_reported_family_change_candidate=int(delta != 0 and reported_bh_pass),
                BH_floor_sensitivity_family_change_candidate=(int(delta != 0 and sensitivity_bh_pass)
                    if sensitivity_bh_pass is not None else None),
                Family_test_scope="family_wide_not_branch_specific"
            ))
        focal_rows.extend(records)
        summary = dict(Species=species, Tested_families=len(tested),
                       Family_p_threshold=args.family_p_threshold,
                       BH_q_threshold=args.fdr_threshold,
                       BH_pvalue_floor_sensitivity=args.pvalue_floor)
        for direction, sign in (("expanded", 1), ("contracted", -1)):
            selected = [row for row in records if row["Copy_change"] * sign > 0]
            summary[f"All_{direction}"] = len(selected)
            summary[f"Nominal_family_p_{direction}"] = sum(row["Family_p_pass"] for row in selected)
            summary[f"BH_reported_{direction}"] = sum(row["Family_BH_reported_pass"] for row in selected)
            summary[f"BH_floor_sensitivity_{direction}"] = (sum(
                row["Family_BH_floor_sensitivity_pass"] or 0 for row in selected)
                if args.pvalue_floor is not None else None)
        summary_rows.append(summary)

    omitted = []
    if args.families:
        with args.families.open(newline="") as handle:
            original = list(csv.DictReader(handle, delimiter="\t"))
        input_ids = {row["Family ID"] for row in original}
        if len(input_ids) != len(original):
            raise ValueError("Duplicate IDs in supplied CAFE family input")
        if set(families) - input_ids:
            raise ValueError("CAFE outputs contain families absent from the supplied input")
        omitted = [dict(Orthogroup=family, Status="Not_in_CAFE_model_output")
                   for family in sorted(input_ids - set(families))]

    args.out.mkdir(parents=True)
    columns = list(focal_rows[0])
    write_table(args.out / "focal_family_changes.tsv", columns, focal_rows)
    write_table(args.out / "focal_nominal_candidates.tsv", columns,
                [row for row in focal_rows if row["Nominal_family_change_candidate"]])
    write_table(args.out / "focal_BH_reported_candidates.tsv", columns,
                [row for row in focal_rows if row["BH_reported_family_change_candidate"]])
    write_table(args.out / "focal_BH_floor_sensitivity_candidates.tsv", columns,
                [row for row in focal_rows if row["BH_floor_sensitivity_family_change_candidate"]])
    write_table(args.out / "focal_summary.tsv", list(summary_rows[0]), summary_rows)
    write_table(args.out / "family_tests.tsv",
                ["Orthogroup", "Family_pvalue", "Family_qvalue_BH_reported",
                 "Family_pvalue_floor_sensitivity", "Family_qvalue_BH_floor_sensitivity"],
                [dict(Orthogroup=family, Family_pvalue=pvalues[family],
                      Family_qvalue_BH_reported=reported_qvalues.get(family),
                      Family_pvalue_floor_sensitivity=sensitivity_pvalues.get(family),
                      Family_qvalue_BH_floor_sensitivity=sensitivity_qvalues.get(family)) for family in ordered])
    write_table(args.out / "family_pvalue_diagnostics.tsv", ["Metric", "Value"], [
        dict(Metric="Tested_families", Value=len(tested)),
        dict(Metric="Reported_zero_pvalues", Value=zero_pvalues),
        dict(Metric="Minimum_positive_reported_pvalue",
             Value=min(positive_pvalues) if positive_pvalues else None),
        dict(Metric="Sensitivity_pvalue_floor", Value=args.pvalue_floor),
        dict(Metric="Families_raised_by_floor", Value=(sum(
            sensitivity_pvalues[family] > value for family, value in tested.items())
            if args.pvalue_floor is not None else None)),
        dict(Metric="Nominal_family_p_threshold", Value=args.family_p_threshold),
        dict(Metric="Supplementary_BH_q_threshold", Value=args.fdr_threshold),
    ])
    write_table(args.out / "input_families_not_in_output.tsv", ["Orthogroup", "Status"], omitted)
    for species in args.species:
        for direction in ("Expansion", "Contraction"):
            for label, flag in [("nominal", "Nominal_family_change_candidate"),
                                ("BH_reported", "BH_reported_family_change_candidate"),
                                ("BH_floor_sensitivity", "BH_floor_sensitivity_family_change_candidate")]:
                identifiers = [row["Orthogroup"] for row in focal_rows
                               if row["Species"] == species and row["Direction"] == direction and row[flag]]
                (args.out / f"{species}.{direction.lower()}.{label}.ids").write_text(
                    "".join(f"{family}\n" for family in identifiers))
    notes = [
        "Schema v2 (2026-11-23): nominal family-wide screening plus reconstructed focal copy changes.",
        "Exploratory results; retain annotation quality, orthogroup definition and time-tree calibration limitations.",
        "BH correction covers all finite family p-values in this model output.",
        "A family-level test is not specific to the focal branch.",
        "A candidate means family-wide p below the chosen threshold AND nonzero reconstructed focal copy change. It is not a branch-specific significance result.",
        "Branch_probability_reported is descriptive only: it is never used as a p-value, threshold or candidate filter. Missing branch probabilities do not exclude a nominal candidate.",
        "For the Gamma model these branch values are posterior probabilities, not p-values: https://github.com/hahnlab/CAFE5/discussions/232",
        "No branch-specific p-values or branch-level FDR conclusions are generated by this script.",
        "Legacy Branch_pvalue_raw / branch_p005 columns and familyFDR_branchP lists are retired; do not interpret old lists as branch-significant results.",
        f"Family p-value floor for BH: {args.pvalue_floor if args.pvalue_floor is not None else 'none'}",
        "Family_pvalue and Family_qvalue_BH_reported preserve reported family p-values and their BH adjustment; any q=0 from reported p=0 remains provisional.",
        "Floor-sensitivity columns are separate and are NA when --pvalue-floor is not supplied. A sensitivity pass is not proof of formal FDR significance.",
        "A p-value floor is a sensitivity assumption; it is not a recovered true p-value or a rigorous confidence bound.",
        "Reported zero family p-values must not be interpreted as exact zero probability; assess simulation and output precision.",
        "If zeros are present, these FDR-based lists are provisional screening candidates pending p-value-resolution review.",
        "Families absent from the model output are not included in the BH denominator.",
        f"Nominal family p threshold: {args.family_p_threshold}",
        f"Supplementary BH q threshold: {args.fdr_threshold}",
        "State how/when the screening threshold was selected; do not call a post-hoc threshold prespecified.",
        f"Results: {args.results.resolve()}",
        f"Tested families: {len(tested)}; untested output families: {len(families) - len(tested)}",
        f"Input families omitted internally: {len(omitted)}",
    ]
    (args.out / "interpretation_notes.txt").write_text("\n".join(notes) + "\n")
    provenance = {kind: dict(path=str(path.resolve()),
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                  for kind, path in paths.items()}
    if args.families:
        provenance["family_input"] = dict(path=str(args.families.resolve()),
            sha256=hashlib.sha256(args.families.read_bytes()).hexdigest())
    (args.out / "input_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (args.out / "summary_settings.json").write_text(json.dumps(dict(
        schema_version=2, family_p_threshold=args.family_p_threshold,
        bh_q_threshold=args.fdr_threshold, bh_pvalue_floor_sensitivity=args.pvalue_floor,
        branch_probabilities_used_as_pvalues=False, model=args.model,
        species=args.species, tested_families=len(tested), reported_zero_pvalues=zero_pvalues,
        raw_results_modified=False), indent=2) + "\n")
    print(f"PASS: {args.out.resolve()}")
    print(f"Family-level tests: {len(tested)}; families absent from CAFE output: {len(omitted)}")
    if zero_pvalues:
        print(f"WARNING: {zero_pvalues} reported family p-values are zero; BH q=0 is not exact-zero evidence. "
              "Review empirical p-value resolution before making formal FDR claims.")
    if args.pvalue_floor is not None:
        print(f"SENSITIVITY ONLY: family p-values floored at {args.pvalue_floor} for BH; "
              "reported p-values are retained unchanged.")
    for row in summary_rows:
        print(row["Species"],
              f"all_expanded={row['All_expanded']}",
              f"all_contracted={row['All_contracted']}",
              f"nominal_family_p_expanded={row['Nominal_family_p_expanded']}",
              f"nominal_family_p_contracted={row['Nominal_family_p_contracted']}",
              f"family_p_threshold={args.family_p_threshold}", sep="\t")


if __name__ == "__main__":
    main()
