#!/usr/bin/env python3
"""Map reviewed protein-ID lists to OrthoFinder orthogroups for exclusion."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_assignment(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Expected Species=/path/to/protein_ids.txt")
    species, path = value.split("=", 1)
    if not species or not path:
        raise argparse.ArgumentTypeError("Expected Species=/path/to/protein_ids.txt")
    return species, Path(path)


def read_ids(path: Path) -> set[str]:
    with path.open() as handle:
        return {
            line.strip().split()[0]
            for line in handle
            if line.strip() and not line.startswith("#")
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", type=Path, required=True,
                        help="OrthoFinder Orthogroups.tsv membership table")
    parser.add_argument("--unassigned", type=Path,
                        help="Orthogroups_UnassignedGenes.tsv; auto-detected when omitted")
    parser.add_argument("--ids", type=parse_assignment, action="append", required=True,
                        help="Repeatable Species=/path/to/protein_ids.txt")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--allow-missing", action="store_true")
    args = parser.parse_args()

    if args.out.exists():
        raise FileExistsError(f"Output already exists: {args.out}")
    args.out.mkdir(parents=True)

    wanted: dict[str, set[str]] = {}
    for species, path in args.ids:
        if species in wanted:
            raise ValueError(f"Species specified more than once: {species}")
        wanted[species] = read_ids(path)
        if not wanted[species]:
            raise ValueError(f"Empty ID list: {path}")

    mappings: list[tuple[str, str, str]] = []
    found = {species: set() for species in wanted}
    family_ids: set[str] = set()

    with args.orthogroups.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup":
            raise ValueError("Unexpected Orthogroups.tsv header")
        absent = set(wanted) - set(reader.fieldnames)
        if absent:
            raise ValueError(f"Species absent from Orthogroups.tsv: {sorted(absent)}")

        for row in reader:
            family = row["Orthogroup"]
            for species, ids in wanted.items():
                members = {x.strip() for x in row[species].split(",") if x.strip()}
                for protein in sorted(ids & members):
                    mappings.append((species, protein, family))
                    found[species].add(protein)
                    family_ids.add(family)

    unresolved = {
        species: sorted(ids - found[species])
        for species, ids in wanted.items()
        if ids - found[species]
    }
    unassigned_path = args.unassigned
    if unassigned_path is None:
        candidate = args.orthogroups.parent / "Orthogroups_UnassignedGenes.tsv"
        if candidate.is_file():
            unassigned_path = candidate

    unassigned_found = {species: set() for species in wanted}
    if unassigned_path is not None:
        if not unassigned_path.is_file():
            raise FileNotFoundError(unassigned_path)
        with unassigned_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup":
                raise ValueError("Unexpected Orthogroups_UnassignedGenes.tsv header")
            absent = set(wanted) - set(reader.fieldnames)
            if absent:
                raise ValueError(
                    f"Species absent from unassigned table: {sorted(absent)}"
                )
            for row in reader:
                for species, ids in wanted.items():
                    members = {
                        x.strip() for x in row[species].split(",") if x.strip()
                    }
                    unassigned_found[species].update(ids & members)

    missing = {
        species: sorted(set(ids) - unassigned_found[species])
        for species, ids in unresolved.items()
        if set(ids) - unassigned_found[species]
    }
    with (args.out / "flagged_protein_to_orthogroup.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Species", "Protein_ID", "Orthogroup"])
        writer.writerows(mappings)
    (args.out / "exclude_orthogroups.ids").write_text(
        "\n".join(sorted(family_ids)) + ("\n" if family_ids else "")
    )
    with (args.out / "unassigned_flagged_proteins.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Species", "Protein_ID", "Interpretation"])
        for species, ids in sorted(unassigned_found.items()):
            writer.writerows(
                (species, protein, "OrthoFinder_unassigned_no_CAFE_family_count")
                for protein in sorted(ids)
            )
    with (args.out / "missing_protein_ids.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Species", "Protein_ID"])
        for species, ids in sorted(missing.items()):
            writer.writerows((species, protein) for protein in ids)

    print(f"Flagged proteins mapped: {len(mappings)}")
    print(f"Orthogroups flagged: {len(family_ids)}")
    print(f"Flagged proteins unassigned by OrthoFinder: "
          f"{sum(map(len, unassigned_found.values()))}")
    print(f"Protein IDs absent from both assigned and unassigned tables: "
          f"{sum(map(len, missing.values()))}")
    if missing and not args.allow_missing:
        raise ValueError("Some protein IDs were not found; inspect missing_protein_ids.tsv")


if __name__ == "__main__":
    main()
