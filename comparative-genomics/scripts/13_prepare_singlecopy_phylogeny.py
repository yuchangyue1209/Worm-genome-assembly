#!/usr/bin/env python3
"""Build a trimmed supermatrix from strict single-copy OrthoFinder groups.

Strict groups are derived from the final Orthogroups.tsv membership table, not
from OrthoFinder's cached SingleCopyOrthologues list. This avoids carrying a
group whose final membership no longer contains exactly one protein per taxon.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            current = raw[1:].split()[0]
            if current in records:
                raise ValueError(f"Duplicate FASTA ID in {path}: {current}")
            records[current] = []
        elif raw.strip():
            if current is None:
                raise ValueError(f"Sequence before header in {path}")
            records[current].append(raw.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def find_sequence(records: dict[str, str], species: str, protein: str) -> str:
    for key in (protein, f"{species}_{protein}"):
        if key in records:
            return records[key]
    suffix = f"_{protein}"
    matches = [seq for name, seq in records.items() if name.endswith(suffix)]
    if len(matches) == 1:
        return matches[0]
    raise KeyError(f"Could not map {species}:{protein} to alignment")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True,
                        help="OrthoFinder Results_* directory")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected-species", type=int, default=11)
    parser.add_argument("--expected-groups", type=int, default=373)
    parser.add_argument("--minimum-column-occupancy", type=int, default=6)
    parser.add_argument("--minimum-columns", type=int, default=50)
    parser.add_argument("--minimum-residues-per-species", type=int, default=20)
    args = parser.parse_args()

    membership = args.results / "Orthogroups/Orthogroups.tsv"
    alignments = args.results / "MultipleSequenceAlignments"
    if not membership.is_file() or not alignments.is_dir():
        raise FileNotFoundError("Missing Orthogroups.tsv or MultipleSequenceAlignments")
    args.out.mkdir(parents=True, exist_ok=False)
    trimmed = args.out / "trimmed"
    trimmed.mkdir()

    with membership.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        species = reader.fieldnames[1:]
        if len(species) != args.expected_species:
            raise ValueError(f"Expected {args.expected_species} species, found {len(species)}")
        groups = []
        for row in reader:
            members = {sp: [x.strip() for x in row[sp].split(",") if x.strip()]
                       for sp in species}
            if all(len(ids) == 1 for ids in members.values()):
                groups.append((row["Orthogroup"], {sp: ids[0] for sp, ids in members.items()}))

    if len(groups) != args.expected_groups:
        raise ValueError(f"Expected {args.expected_groups} strict groups, found {len(groups)}")

    supermatrix = {sp: [] for sp in species}
    partitions = []
    audit = []
    position = 1
    retained_ids = []
    for orthogroup, members in groups:
        source = alignments / f"{orthogroup}.fa"
        records = read_fasta(source)
        seqs = {sp: find_sequence(records, sp, pid) for sp, pid in members.items()}
        lengths = {len(seq) for seq in seqs.values()}
        if len(lengths) != 1:
            raise ValueError(f"Unequal alignment lengths for {orthogroup}: {sorted(lengths)}")
        length = lengths.pop()
        keep = [i for i in range(length)
                if sum(seqs[sp][i] in STANDARD for sp in species) >= args.minimum_column_occupancy]
        filtered = {sp: "".join(seqs[sp][i] if seqs[sp][i] in STANDARD else "X" for i in keep)
                    for sp in species}
        residues = {sp: sum(aa in STANDARD for aa in seq) for sp, seq in filtered.items()}
        accepted = (len(keep) >= args.minimum_columns and
                    min(residues.values(), default=0) >= args.minimum_residues_per_species)
        audit.append((orthogroup, length, len(keep), min(residues.values()), int(accepted)))
        if not accepted:
            continue
        retained_ids.append(orthogroup)
        (trimmed / f"{orthogroup}.fa").write_text(
            "".join(f">{sp}\n{filtered[sp]}\n" for sp in species))
        for sp in species:
            supermatrix[sp].append(filtered[sp])
        end = position + len(keep) - 1
        partitions.append((orthogroup, position, end))
        position = end + 1

    (args.out / "supermatrix.fa").write_text(
        "".join(f">{sp}\n{''.join(supermatrix[sp])}\n" for sp in species))
    (args.out / "partitions.nex").write_text(
        "#nexus\nbegin sets;\n" +
        "".join(f"  charset {og} = {start}-{end};\n" for og, start, end in partitions) +
        "end;\n")
    (args.out / "strict_singlecopy_from_final_membership.txt").write_text(
        "\n".join(og for og, _ in groups) + "\n")
    (args.out / "retained_loci.txt").write_text("\n".join(retained_ids) + "\n")
    with (args.out / "locus_filter_audit.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Orthogroup", "Original_columns", "Retained_columns",
                         "Minimum_species_residues", "Accepted"])
        writer.writerows(audit)
    with (args.out / "sha256.tsv").open("w") as handle:
        for name in ("supermatrix.fa", "partitions.nex", "retained_loci.txt"):
            path = args.out / name
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}\t{name}\n")
    print(f"Strict final-membership groups: {len(groups)}")
    print(f"Retained loci: {len(retained_ids)}")
    print(f"Supermatrix sites: {position - 1}")


if __name__ == "__main__":
    main()
