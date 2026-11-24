#!/usr/bin/env python3
"""Integrate focal species-specific orthogroups with reciprocal homology evidence.

The output distinguishes raw OrthoFinder species-specific orthogroups from
proteins that remain unresolved after sister-proteome BLASTP and sister-genome
TBLASTN searches.  InterPro and reviewed TE-candidate lists are optional QC
layers; they do not redefine OrthoFinder membership.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FOCAL_PAIRS = {
    "S_cotti": "S_solidus",
    "S_solidus": "S_cotti",
}


def key_value(text: str) -> tuple[str, Path]:
    key, value = text.split("=", 1)
    return key, Path(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--search-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--qcov", type=float, default=50.0)
    parser.add_argument(
        "--interpro",
        action="append",
        default=[],
        metavar="SPECIES=TSV",
    )
    parser.add_argument(
        "--te-ids",
        action="append",
        default=[],
        metavar="SPECIES=IDS",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    parts: list[str] = []

    def save() -> None:
        if header is None:
            return
        protein_id = header.split()[0]
        sequence = "".join(parts)
        if protein_id in records or not sequence:
            raise ValueError(f"Invalid FASTA record in {path}: {protein_id}")
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
                    raise ValueError(f"Sequence before header in {path}")
                parts.append(line)
    save()
    return records


def best_qcov(path: Path, qcov_column: int) -> dict[str, float]:
    result = {}
    if not path.is_file() or not path.stat().st_size:
        return result
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            result[fields[0]] = float(fields[qcov_column])
    return result


def read_annotations(
    specifications: list[str], wanted: set[tuple[str, str]]
) -> dict[tuple[str, str], dict[str, set[str]]]:
    annotations: dict[tuple[str, str], dict[str, set[str]]] = {}
    for specification in specifications:
        species, path = key_value(specification)
        if not path.is_file():
            raise FileNotFoundError(path)
        wanted_ids = {protein for sp, protein in wanted if sp == species}
        with path.open(errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[0] not in wanted_ids:
                    continue
                record = annotations.setdefault(
                    (species, fields[0]),
                    {"pfam": set(), "interpro": set(), "descriptions": set()},
                )
                if fields[3] == "Pfam":
                    record["pfam"].add(fields[4])
                if len(fields) > 11 and fields[11] not in {"", "-", "NULL"}:
                    record["interpro"].add(fields[11])
                for index in (5, 12):
                    if len(fields) > index and fields[index] not in {"", "-", "NULL"}:
                        record["descriptions"].add(fields[index])
    return annotations


def write_subset(
    output: Path, protein_ids: list[str], records: dict[str, tuple[str, str]]
) -> None:
    with output.open("x") as handle:
        for protein_id in protein_ids:
            header, sequence = records[protein_id]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    args = parse_args()
    if args.out.exists():
        raise FileExistsError(
            f"Output directory already exists: {args.out}. Use a versioned path."
        )
    args.out.mkdir(parents=True)

    membership_path = args.candidates / "focal_species_specific_orthogroups.tsv"
    membership: dict[tuple[str, str], str] = {}
    members_by_group: dict[tuple[str, str], list[str]] = {}
    with membership_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            proteins = [item.strip() for item in row["Protein_IDs"].split(",") if item.strip()]
            members_by_group[(row["Species"], row["Orthogroup"])] = proteins
            for protein in proteins:
                membership[(row["Species"], protein)] = row["Orthogroup"]

    fasta = {
        species: read_fasta(args.candidates / f"{species}.species_specific.faa")
        for species in FOCAL_PAIRS
    }
    all_proteins = {
        (species, protein)
        for species, records in fasta.items()
        for protein in records
    }

    te_ids: dict[str, set[str]] = defaultdict(set)
    for specification in args.te_ids:
        species, path = key_value(specification)
        if not path.is_file():
            raise FileNotFoundError(path)
        te_ids[species] = {line.strip() for line in path.open() if line.strip()}

    evidence_rows = []
    unresolved: dict[str, list[str]] = defaultdict(list)
    group_counts: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: {"protein": 0, "genome": 0, "unresolved": 0}
    )

    for species, target in FOCAL_PAIRS.items():
        protein_hits = best_qcov(
            args.search_dir / "blastp" / f"{species}_vs_{target}.best.tsv", 6
        )
        genome_hits = best_qcov(
            args.search_dir / "tblastn" / f"{species}_vs_{target}_genome.best.tsv", 9
        )
        for protein_id, (_, sequence) in fasta[species].items():
            orthogroup = membership[(species, protein_id)]
            protein_qcov = protein_hits.get(protein_id, 0.0)
            genome_tested = protein_qcov < args.qcov
            genome_qcov = genome_hits.get(protein_id, 0.0) if genome_tested else 0.0
            if protein_qcov >= args.qcov:
                evidence = "Sister_proteome"
                group_counts[(species, orthogroup)]["protein"] += 1
            elif genome_qcov >= args.qcov:
                evidence = "Sister_genome"
                group_counts[(species, orthogroup)]["genome"] += 1
            else:
                evidence = "Unresolved"
                unresolved[species].append(protein_id)
                group_counts[(species, orthogroup)]["unresolved"] += 1
            evidence_rows.append(
                {
                    "Species": species,
                    "Orthogroup": orthogroup,
                    "Protein_ID": protein_id,
                    "Length_aa": len(sequence),
                    "Protein_qcov": f"{protein_qcov:.1f}",
                    "Genome_tested": int(genome_tested),
                    "Genome_qcov": f"{genome_qcov:.1f}",
                    "Cross_species_evidence": evidence,
                }
            )

    with (args.out / "species_specific_final_evidence.tsv").open(
        "x", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(evidence_rows)

    with (args.out / "species_specific_orthogroup_evidence.tsv").open(
        "x", newline=""
    ) as handle:
        fields = [
            "Species",
            "Orthogroup",
            "Proteins",
            "Sister_proteome_evidence",
            "Sister_genome_evidence",
            "Unresolved",
            "Fully_unresolved_orthogroup",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (species, orthogroup), proteins in members_by_group.items():
            counts = group_counts[(species, orthogroup)]
            writer.writerow(
                {
                    "Species": species,
                    "Orthogroup": orthogroup,
                    "Proteins": len(proteins),
                    "Sister_proteome_evidence": counts["protein"],
                    "Sister_genome_evidence": counts["genome"],
                    "Unresolved": counts["unresolved"],
                    "Fully_unresolved_orthogroup": int(
                        counts["unresolved"] == len(proteins)
                    ),
                }
            )

    for species in FOCAL_PAIRS:
        ids = unresolved[species]
        (args.out / f"{species}.unresolved_after_protein_and_genome.ids").write_text(
            "".join(f"{protein_id}\n" for protein_id in ids)
        )
        write_subset(
            args.out / f"{species}.unresolved_after_protein_and_genome.faa",
            ids,
            fasta[species],
        )

    unresolved_set = {
        (species, protein) for species, proteins in unresolved.items() for protein in proteins
    }
    annotations = read_annotations(args.interpro, unresolved_set)
    review_rows = []
    for species, proteins in unresolved.items():
        for protein_id in proteins:
            sequence = fasta[species][protein_id][1]
            annotation = annotations.get(
                (species, protein_id),
                {"pfam": set(), "interpro": set(), "descriptions": set()},
            )
            flags = []
            if protein_id in te_ids[species]:
                flags.append("STRONG_TE")
            if len(sequence) < 50:
                flags.append("VERY_SHORT_LT50")
            elif len(sequence) < 100:
                flags.append("SHORT_LT100")
            if annotation["pfam"] or annotation["interpro"]:
                flags.append("DOMAIN_SUPPORTED")
            else:
                flags.append("NO_DOMAIN")
            if "STRONG_TE" in flags:
                priority = "Exclude_TE"
            elif "VERY_SHORT_LT50" in flags:
                priority = "Low"
            elif "DOMAIN_SUPPORTED" in flags:
                priority = "High"
            else:
                priority = "Medium"
            review_rows.append(
                {
                    "Species": species,
                    "Orthogroup": membership[(species, protein_id)],
                    "Protein_ID": protein_id,
                    "Length_aa": len(sequence),
                    "Pfam": ";".join(sorted(annotation["pfam"])) or "NA",
                    "InterPro": ";".join(sorted(annotation["interpro"])) or "NA",
                    "Descriptions": ";".join(sorted(annotation["descriptions"])) or "NA",
                    "Flags": ";".join(flags),
                    "Review_priority": priority,
                }
            )

    with (args.out / "unresolved_candidate_review.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(review_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(review_rows)

    for species in FOCAL_PAIRS:
        groups = {
            membership[(species, protein_id)] for protein_id in unresolved[species]
        }
        print(
            f"{species}\tunresolved_proteins={len(unresolved[species])}"
            f"\torthogroups_with_unresolved_members={len(groups)}"
        )
    print(f"PASS: {args.out}")


if __name__ == "__main__":
    main()
