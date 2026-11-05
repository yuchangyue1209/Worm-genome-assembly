#!/usr/bin/env python3
"""Audit SprJ_4.0 gene, transcript, CDS and protein identifier relationships."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote


def parse_attributes(text: str) -> dict[str, str]:
    return dict(field.split("=", 1) for field in text.split(";") if "=" in field)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.gff, args.proteins):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing or empty {path}")
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: output directory already exists: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    proteins: set[str] = set()
    for line in args.proteins.open():
        if line.startswith(">"):
            protein_id = line[1:].split()[0]
            if protein_id in proteins:
                raise SystemExit(f"ERROR: duplicate protein ID: {protein_id}")
            proteins.add(protein_id)

    genes: set[str] = set()
    mrnas: dict[str, tuple[str, ...]] = {}
    cds_count: Counter[str] = Counter()
    protein_ids: defaultdict[str, set[str]] = defaultdict(set)
    transcripts_per_gene: Counter[str] = Counter()

    for line in args.gff.open():
        if line.startswith("##FASTA"):
            break
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) != 9:
            continue
        attrs = parse_attributes(fields[8])
        if fields[2] == "gene":
            genes.add(attrs["ID"])
        elif fields[2] in {"mRNA", "transcript"}:
            transcript_id = attrs["ID"]
            if transcript_id in mrnas:
                raise SystemExit(f"ERROR: duplicate transcript ID: {transcript_id}")
            parent = attrs.get("Parent", "")
            mrnas[transcript_id] = (
                parent,
                fields[0],
                fields[3],
                fields[4],
                fields[6],
                fields[8],
            )
            transcripts_per_gene[parent] += 1
        elif fields[2] == "CDS":
            for parent in attrs.get("Parent", "").split(","):
                if parent:
                    cds_count[parent] += 1
                    if attrs.get("protein_id"):
                        protein_ids[parent].add(unquote(attrs["protein_id"]))

    status_counts: Counter[str] = Counter()
    output = args.output_dir / "transcript_protein_audit.tsv"
    with output.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "transcript_id",
                "gene_id",
                "seqid",
                "start",
                "end",
                "strand",
                "CDS_features",
                "protein_ids",
                "status",
                "attributes",
            ]
        )
        for transcript_id, record in sorted(mrnas.items()):
            ids = protein_ids[transcript_id]
            if not cds_count[transcript_id]:
                status = "no_CDS"
            elif not ids:
                status = "CDS_without_protein_id"
            elif not ids <= proteins:
                status = "protein_id_missing_from_FASTA"
            else:
                status = "protein_present"
            status_counts[status] += 1
            writer.writerow(
                [
                    transcript_id,
                    *record[:5],
                    cds_count[transcript_id],
                    ",".join(sorted(ids)),
                    status,
                    record[5],
                ]
            )

    referenced = set().union(*protein_ids.values()) if protein_ids else set()
    unreferenced = proteins - referenced
    (args.output_dir / "proteins_not_referenced_by_CDS.txt").write_text(
        "".join(f"{protein_id}\n" for protein_id in sorted(unreferenced))
    )

    print(f"Genes: {len(genes)}")
    print(f"Transcripts: {len(mrnas)}")
    print(f"Protein FASTA sequences: {len(proteins)}")
    print(
        "Genes with multiple transcripts:",
        sum(count > 1 for count in transcripts_per_gene.values()),
    )
    for status, count in sorted(status_counts.items()):
        print(f"{status}: {count}")
    print(f"Proteins not referenced by CDS: {len(unreferenced)}")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
