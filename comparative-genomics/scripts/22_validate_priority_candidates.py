#!/usr/bin/env python3
"""Validate manually prioritized focal candidates with RNA and repeat evidence.

The candidate list is derived from an explicit orthogroup review table.  CDS
coordinates are extracted from each structural annotation, RNA depth is
calculated with samtools, and CDS overlap with a RepeatMasker GFF is reported
when available.  RNA support validates transcription, not species specificity.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--priority-table", required=True, type=Path)
    parser.add_argument("--datasets", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--mapq", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    attributes = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def load_priority(path: Path) -> set[tuple[str, str]]:
    selected = set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Species", "Orthogroup", "Use_for_RNA_repeat_validation"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Priority table requires columns: {sorted(required)}")
        for row in reader:
            if row["Use_for_RNA_repeat_validation"] == "1":
                selected.add((row["Species"], row["Orthogroup"]))
    if not selected:
        raise ValueError("Priority table selected no orthogroups")
    return selected


def load_members(
    path: Path, selected: set[tuple[str, str]]
) -> dict[str, dict[str, str]]:
    members: dict[str, dict[str, str]] = defaultdict(dict)
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["Species"], row["Orthogroup"])
            if key not in selected:
                continue
            for protein in row["Protein_IDs"].split(","):
                members[row["Species"]][protein.strip()] = row["Orthogroup"]
    missing = selected - {
        (species, orthogroup)
        for species, mapping in members.items()
        for orthogroup in mapping.values()
    }
    if missing:
        raise ValueError(f"Priority orthogroups absent from membership table: {missing}")
    return members


def load_datasets(path: Path) -> dict[str, dict[str, Path | None]]:
    datasets = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Species", "GFF3", "RNA_BAM", "Repeat_GFF"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Dataset table requires columns: {sorted(required)}")
        for row in reader:
            repeat = None if row["Repeat_GFF"] in {"", "NA", "-"} else Path(row["Repeat_GFF"])
            datasets[row["Species"]] = {
                "gff": Path(row["GFF3"]),
                "bam": Path(row["RNA_BAM"]),
                "repeat": repeat,
            }
    return datasets


def extract_cds(
    gff: Path, wanted: set[str]
) -> dict[str, list[tuple[str, int, int, str]]]:
    result: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    with gff.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_attributes(fields[8])
            parents = set(attrs.get("Parent", "").split(","))
            protein_id = attrs.get("protein_id", "")
            for candidate in wanted:
                if candidate in parents or candidate == protein_id:
                    result[candidate].append(
                        (fields[0], int(fields[3]), int(fields[4]), fields[6])
                    )
    missing = sorted(wanted - set(result))
    if missing:
        raise ValueError(f"CDS not found in {gff}: {missing}")
    return result


def merged_cds(
    raw: dict[str, list[tuple[str, int, int, str]]]
) -> dict[str, list[tuple[str, int, int, str]]]:
    result = {}
    for protein, records in raw.items():
        groups: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
        for chrom, start, end, strand in records:
            groups[(chrom, strand)].append((start, end))
        merged = []
        for (chrom, strand), intervals in groups.items():
            for start, end in merge_intervals(intervals):
                merged.append((chrom, start, end, strand))
        result[protein] = merged
    return result


def repeat_overlap(
    repeat_gff: Path | None,
    cds: dict[str, list[tuple[str, int, int, str]]],
) -> dict[str, int | None]:
    if repeat_gff is None:
        return {protein: None for protein in cds}
    if not repeat_gff.is_file():
        raise FileNotFoundError(repeat_gff)
    by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for protein, intervals in cds.items():
        for chrom, start, end, _ in intervals:
            by_chrom[chrom].append((start, end, protein))
    overlap_segments: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with repeat_gff.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 5 or fields[0] not in by_chrom:
                continue
            repeat_start, repeat_end = int(fields[3]), int(fields[4])
            for start, end, protein in by_chrom[fields[0]]:
                left = max(start, repeat_start)
                right = min(end, repeat_end)
                if left <= right:
                    overlap_segments[protein].append((left, right))
    result: dict[str, int | None] = {}
    for protein in cds:
        result[protein] = sum(
            end - start + 1
            for start, end in merge_intervals(overlap_segments[protein])
        )
    return result


def main() -> None:
    args = parse_args()
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools is required")
    if args.out.exists() and any(args.out.iterdir()) and not args.force:
        raise FileExistsError(f"Non-empty output exists: {args.out}; use --force or a new path")
    args.out.mkdir(parents=True, exist_ok=True)

    selected = load_priority(args.priority_table)
    members = load_members(args.membership, selected)
    datasets = load_datasets(args.datasets)
    output_rows = []

    for species, protein_to_og in members.items():
        if species not in datasets:
            raise ValueError(f"No validation dataset for {species}")
        dataset = datasets[species]
        gff = dataset["gff"]
        bam = dataset["bam"]
        if not isinstance(gff, Path) or not gff.is_file():
            raise FileNotFoundError(gff)
        if not isinstance(bam, Path) or not bam.is_file():
            raise FileNotFoundError(bam)
        subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)

        cds = merged_cds(extract_cds(gff, set(protein_to_og)))
        bed = args.out / f"{species}.priority_candidate_cds.bed"
        with bed.open("w") as handle:
            for protein in sorted(cds):
                for chrom, start, end, strand in cds[protein]:
                    handle.write(
                        f"{chrom}\t{start - 1}\t{end}\t{protein}\t0\t{strand}\n"
                    )

        depth_file = args.out / f"{species}.priority_candidate.depth.tsv"
        with depth_file.open("w") as handle:
            subprocess.run(
                [
                    "samtools",
                    "depth",
                    "-a",
                    "-Q",
                    str(args.mapq),
                    "-b",
                    str(bed),
                    str(bam),
                ],
                check=True,
                stdout=handle,
            )

        intervals_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        cds_bp = {}
        for protein, intervals in cds.items():
            cds_bp[protein] = sum(end - start + 1 for _, start, end, _ in intervals)
            for chrom, start, end, _ in intervals:
                intervals_by_chrom[chrom].append((start, end, protein))

        depth_sum = defaultdict(int)
        covered1 = defaultdict(int)
        covered5 = defaultdict(int)
        with depth_file.open() as handle:
            for line in handle:
                chrom, position_text, depth_text = line.rstrip().split("\t")[:3]
                position, depth = int(position_text), int(depth_text)
                for start, end, protein in intervals_by_chrom.get(chrom, []):
                    if start <= position <= end:
                        depth_sum[protein] += depth
                        covered1[protein] += depth >= 1
                        covered5[protein] += depth >= 5

        overlap = repeat_overlap(dataset["repeat"], cds)
        for protein in sorted(protein_to_og):
            length = cds_bp[protein]
            mean_depth = depth_sum[protein] / length
            pct1 = 100 * covered1[protein] / length
            pct5 = 100 * covered5[protein] / length
            if pct1 >= 80 and pct5 >= 50:
                support = "Strong"
            elif pct1 >= 50:
                support = "Moderate"
            elif pct1 >= 20:
                support = "Partial"
            else:
                support = "Weak_or_not_detected"
            repeat_bp = overlap[protein]
            output_rows.append(
                {
                    "Species": species,
                    "Orthogroup": protein_to_og[protein],
                    "Protein_ID": protein,
                    "CDS_bp": length,
                    "Mean_depth": f"{mean_depth:.2f}",
                    "Covered_ge1x_percent": f"{pct1:.1f}",
                    "Covered_ge5x_percent": f"{pct5:.1f}",
                    "RNA_support": support,
                    "Repeat_overlap_bp": "NA" if repeat_bp is None else repeat_bp,
                    "Repeat_overlap_percent": (
                        "NA" if repeat_bp is None else f"{100 * repeat_bp / length:.2f}"
                    ),
                }
            )

    output = args.out / "priority_candidate_rnaseq_repeat_support.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"PASS: {output}")
    print(f"Proteins summarized: {len(output_rows)}")


if __name__ == "__main__":
    main()
