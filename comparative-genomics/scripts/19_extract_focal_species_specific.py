#!/usr/bin/env python3
"""Extract focal-species-only OrthoFinder orthogroups and their proteins.

The script reads the final membership table (Orthogroups.tsv), not cached
single-copy or species-specific lists.  A focal orthogroup is considered raw
species-specific when it contains one or more proteins from exactly one
species and no proteins from any other species in that OrthoFinder run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--species",
        nargs="+",
        default=["S_cotti", "S_solidus"],
        help="Focal species names matching Orthogroups.tsv columns",
    )
    parser.add_argument(
        "--expected",
        action="append",
        default=[],
        metavar="SPECIES=COUNT",
        help="Optional raw species-specific orthogroup count assertion",
    )
    return parser.parse_args()


def split_members(text: str) -> list[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    parts: list[str] = []

    def save() -> None:
        if header is None:
            return
        protein_id = header.split()[0]
        if protein_id in records:
            raise ValueError(f"Duplicate FASTA ID in {path}: {protein_id}")
        sequence = "".join(parts)
        if not sequence:
            raise ValueError(f"Empty FASTA sequence in {path}: {protein_id}")
        records[protein_id] = (header, sequence)

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                header = line[1:]
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before first header in {path}")
                parts.append(line)
    save()
    return records


def fasta_path(input_dir: Path, species: str) -> Path:
    candidates = [
        input_dir / f"{species}.faa",
        input_dir / f"{species}.fa",
        input_dir / f"{species}.fasta",
    ]
    found = [path for path in candidates if path.is_file() and path.stat().st_size]
    if len(found) != 1:
        raise FileNotFoundError(
            f"Expected one non-empty FASTA for {species} in {input_dir}; found {found}"
        )
    return found[0]


def write_fasta(path: Path, ids: list[str], records: dict[str, tuple[str, str]]) -> None:
    missing = [protein_id for protein_id in ids if protein_id not in records]
    if missing:
        raise KeyError(f"Missing {len(missing)} proteins in FASTA {path}: {missing[:5]}")
    with path.open("x") as handle:
        for protein_id in ids:
            header, sequence = records[protein_id]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    args = parse_args()
    orthogroups = args.results / "Orthogroups" / "Orthogroups.tsv"
    if not orthogroups.is_file():
        raise FileNotFoundError(orthogroups)
    if args.out.exists():
        raise FileExistsError(
            f"Output directory already exists: {args.out}. Use a new versioned directory."
        )

    expected: dict[str, int] = {}
    for item in args.expected:
        species, count = item.split("=", 1)
        expected[species] = int(count)

    args.out.mkdir(parents=True)
    rows: list[dict[str, str | int]] = []
    ids_by_species: dict[str, list[str]] = defaultdict(list)

    with orthogroups.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup":
            raise ValueError(f"Unexpected Orthogroups.tsv header: {reader.fieldnames}")
        all_species = reader.fieldnames[1:]
        missing_species = sorted(set(args.species) - set(all_species))
        if missing_species:
            raise ValueError(f"Focal species absent from OrthoFinder results: {missing_species}")

        for row in reader:
            members = {species: split_members(row[species]) for species in all_species}
            present = [species for species, proteins in members.items() if proteins]
            if len(present) != 1 or present[0] not in args.species:
                continue
            species = present[0]
            proteins = members[species]
            rows.append(
                {
                    "Species": species,
                    "Orthogroup": row["Orthogroup"],
                    "Protein_count": len(proteins),
                    "Protein_IDs": ",".join(proteins),
                }
            )
            ids_by_species[species].extend(proteins)

    summary = args.out / "focal_species_specific_orthogroups.tsv"
    with summary.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Species", "Orthogroup", "Protein_count", "Protein_IDs"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    length_rows: list[tuple[str, str, int]] = []
    manifest = []
    for species in args.species:
        species_rows = [row for row in rows if row["Species"] == species]
        if species in expected and len(species_rows) != expected[species]:
            raise ValueError(
                f"{species}: expected {expected[species]} species-specific orthogroups; "
                f"found {len(species_rows)}"
            )
        source = fasta_path(args.input_dir, species)
        records = read_fasta(source)
        ids = ids_by_species[species]
        (args.out / f"{species}.species_specific.ids").write_text(
            "".join(f"{protein_id}\n" for protein_id in ids)
        )
        output_fasta = args.out / f"{species}.species_specific.faa"
        write_fasta(output_fasta, ids, records)
        for protein_id in ids:
            length_rows.append((species, protein_id, len(records[protein_id][1])))
        manifest.append(
            {
                "species": species,
                "orthogroups": len(species_rows),
                "proteins": len(ids),
                "source": str(source.resolve()),
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            }
        )
        print(f"{species}\torthogroups={len(species_rows)}\tproteins={len(ids)}")

    with (args.out / "focal_species_specific_protein_lengths.tsv").open(
        "x", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Species", "Protein_ID", "Length_aa"])
        writer.writerows(length_rows)

    metadata = {
        "results": str(args.results.resolve()),
        "orthogroups_tsv_sha256": hashlib.sha256(orthogroups.read_bytes()).hexdigest(),
        "species": args.species,
        "inputs": manifest,
    }
    (args.out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"PASS: {summary}")


if __name__ == "__main__":
    main()
