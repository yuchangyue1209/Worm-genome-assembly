#!/usr/bin/env python3
"""Build the next Scotti TE-review packet without changing any gene models."""
from pathlib import Path
import argparse
import csv


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--review",
        type=Path,
        default=Path(
            "/work/cyu/annotation/scotti-annotation/results/"
            "annotation_qc/08_TE_protein_review_v1"
        ),
    )
    args = parser.parse_args()
    root = args.review
    first = root / "first30_case_review_v1"

    previous = read_table(first / "first30.primary_case_review.tsv")
    pending = [
        row for row in previous
        if row["case_review_category"] not in {
            "strong_TE_associated_candidate", "cellular_function_supported"
        }
    ]
    remaining = read_table(first / "remaining_functional_candidates.tsv")
    rows = pending + remaining
    assert len(pending) == 14, len(pending)
    assert len(remaining) == 70, len(remaining)
    assert len({row["protein_id"] for row in rows}) == 84

    domains = {}
    for row in read_table(root / "matched_interpro_evidence.tsv"):
        domains.setdefault(row["protein_id"], []).append(row)

    fields = (
        "length_aa", "CDS_RM_classified_TE_pct", "CDS_RM_unknown_pct",
        "explicit_TE_signals", "RT_integrase_signals", "generic_signals",
        "host_context_signals", "swissprot_title", "swissprot_evalue",
        "swissprot_query_coverage_pct", "eggnog_description",
        "interpro_descriptions",
    )
    target = first / "next84_review_packet.txt"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite: {target}")
    with target.open("x") as handle:
        for row in rows:
            handle.write(f"\n=== {row['gene_id']} / {row['protein_id']} ===\n")
            handle.write(
                "review_category: "
                + row.get("case_review_category", row["review_tier"])
                + "\n"
            )
            for key in fields:
                handle.write(f"{key}: {row.get(key) or 'none'}\n")
            handle.write("Matching member-database evidence:\n")
            for domain in domains.get(row["protein_id"], []):
                values = [
                    domain["database"], domain["signature"],
                    domain["signature_description"],
                    f"{domain['aa_start']}-{domain['aa_end']}",
                    f"score={domain['score']}",
                ]
                handle.write("  " + " | ".join(values) + "\n")
    print("Unresolved previous models:", len(pending))
    print("Remaining functional candidates:", len(remaining))
    print("Packet:", target)


if __name__ == "__main__":
    main()
