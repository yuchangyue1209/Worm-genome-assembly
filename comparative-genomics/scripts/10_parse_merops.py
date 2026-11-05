#!/usr/bin/env python3
"""Filter MEROPS BLASTP hits, select best hits, and attach metadata."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

HEADER_RE = re.compile(r"^>\s*(MER\d+)\s+-\s+(.*?)\s+\[([^\]]+)\]#([^#]+)#")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--merops-fasta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--subject-coverage", type=float, default=50.0)
    parser.add_argument("--species", nargs="+", required=True)
    return parser.parse_args()


def metadata(path: Path) -> dict[str, tuple[str, str, str]]:
    result: dict[str, tuple[str, str, str]] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = HEADER_RE.match(line.rstrip())
            if match:
                accession, name, merops_id, family = match.groups()
                result[accession] = (name, merops_id, family)
    if not result:
        raise ValueError(f"No MEROPS metadata parsed from {path}")
    return result


def merops_class(family: str) -> str:
    return "Protease_inhibitor" if family.upper().startswith("I") else "Protease"


def rank(hit: dict[str, object]) -> tuple[float, float, float, float]:
    return (
        float(hit["bitscore"]),
        -float(hit["evalue"]),
        float(hit["subject_coverage"]),
        float(hit["pident"]),
    )


def select_best(
    path: Path, max_evalue: float, min_subject_coverage: float
) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    with path.open() as handle:
        for number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields:
                continue
            if len(fields) != 14:
                raise ValueError(f"{path}:{number}: expected 14 columns, found {len(fields)}")
            query, subject, pident, length = fields[:4]
            evalue, bitscore, slen = fields[10], fields[11], fields[13]
            try:
                # Preserve the legacy project definition recorded in the
                # confirmed table: aligned residues divided by subject length.
                coverage = 100.0 * int(length) / int(slen)
                hit: dict[str, object] = {
                    "subject": subject,
                    "pident": float(pident),
                    "length": int(length),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                    "subject_coverage": coverage,
                }
            except (ValueError, ZeroDivisionError) as error:
                raise ValueError(f"{path}:{number}: malformed BLAST row") from error
            if float(hit["evalue"]) > max_evalue or coverage < min_subject_coverage:
                continue
            current = best.get(query)
            if current is None or rank(hit) > rank(current):
                best[query] = hit
    return best


def fmt(value: float) -> str:
    return f"{value:.6g}"


def main() -> None:
    args = arguments()
    annotations = metadata(args.merops_fasta)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[tuple[str, str]] = Counter()
    rows: list[list[object]] = []

    for species in args.species:
        raw = args.raw_dir / f"{species}.merops.blastp.tsv"
        if not raw.is_file() or raw.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty MEROPS result: {raw}")
        for query, hit in select_best(raw, args.evalue, args.subject_coverage).items():
            accession = str(hit["subject"])
            if accession not in annotations:
                raise KeyError(f"MEROPS accession absent from pepunit.lib: {accession}")
            name, merops_id, family = annotations[accession]
            category = merops_class(family)
            counts[(species, category)] += 1
            rows.append([
                species, query, accession, name, merops_id, family, category,
                fmt(float(hit["pident"])), hit["length"], fmt(float(hit["evalue"])),
                fmt(float(hit["bitscore"])), f"{float(hit['subject_coverage']):.1f}",
            ])

    with args.output.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "Species", "Protein_ID", "MEROPS_accession", "MEROPS_name", "MEROPS_ID",
            "MEROPS_family", "MEROPS_class", "Percent_identity", "Alignment_length",
            "Evalue", "Bitscore", "Subject_coverage",
        ])
        writer.writerows(rows)

    with args.summary.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["Species", "Protease", "Protease_inhibitor", "Total"])
        for species in args.species:
            protease = counts[(species, "Protease")]
            inhibitor = counts[(species, "Protease_inhibitor")]
            writer.writerow([species, protease, inhibitor, protease + inhibitor])

    print(f"PASS: {args.output}")
    print(f"PASS: {args.summary}")


if __name__ == "__main__":
    main()
