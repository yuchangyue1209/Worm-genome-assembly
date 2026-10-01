#!/usr/bin/env python3
"""Select extreme OrthoFinder families and flag likely TE-derived groups."""

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

TE = re.compile(
    r"reverse transcriptase|retrotranspos|retroelement|transposase|transposition|"
    r"integrase|ribonuclease h|rna-directed dna polymerase|mobile element|"
    r"PF00078|PF00665|PF14529|PF03372|PF17921|PF17917",
    re.I,
)


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path,
                        help="OrthoFinder results directory")
    parser.add_argument("--interpro", required=True, type=Path)
    parser.add_argument("--eggnog", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--minimum-copy", type=int, default=50)
    return parser.parse_args()


def main():
    args = arguments()
    args.out.mkdir(parents=True, exist_ok=True)
    counts_file = args.results / "Orthogroups/Orthogroups.GeneCount.tsv"
    groups_file = args.results / "Orthogroups/Orthogroups.tsv"

    with counts_file.open(newline="") as handle:
        rows = list(csv.DictReader((x.replace("\r", "") for x in handle), delimiter="\t"))
    species = [x for x in rows[0] if x not in {"Orthogroup", "Total"}]
    selected = {
        row["Orthogroup"] for row in rows
        if max(int(row[x]) for x in species) >= args.minimum_copy
    }

    members = defaultdict(list)
    wanted = set()
    with groups_file.open(newline="") as handle:
        reader = csv.reader((x.replace("\r", "") for x in handle), delimiter="\t")
        header = next(reader)
        for row in reader:
            if row[0] not in selected:
                continue
            for sp, field in zip(header[1:], row[1:]):
                for protein in filter(None, (x.strip() for x in field.split(", "))):
                    members[row[0]].append((sp, protein))
                    wanted.add((sp, protein))

    annotation = defaultdict(lambda: {"egg": set(), "ipr": set(), "pfam": set()})
    for sp in species:
        with (args.interpro / f"{sp}.tsv").open(errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 5 or (sp, fields[0]) not in wanted:
                    continue
                key = (sp, fields[0])
                if fields[3] == "Pfam": annotation[key]["pfam"].add(fields[4])
                if len(fields) >= 13 and fields[11].startswith("IPR"):
                    annotation[key]["ipr"].add(f"{fields[11]}:{fields[12]}")
        with (args.eggnog / f"{sp}.emapper.annotations").open(errors="replace") as handle:
            for line in handle:
                if line.startswith("#"): continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) >= 8 and (sp, fields[0]) in wanted and fields[7] != "-":
                    annotation[(sp, fields[0])]["egg"].add(fields[7])

    output = args.out / "extreme_family_annotation_summary.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Orthogroup", "Total_members", "Annotated_members",
                         "Annotation_percent", "TE_keyword_hits", "Initial_classification",
                         "Top_eggNOG", "Top_InterPro", "Top_Pfam"])
        for og in sorted(selected):
            counters = {key: Counter() for key in ("egg", "ipr", "pfam")}
            annotated = 0
            for item in members[og]:
                values = annotation[item]
                annotated += int(any(values.values()))
                for key in counters:
                    counters[key].update(values[key])
            labels = " ".join(x for counter in counters.values() for x in counter)
            hits = sorted({x.group(0).lower() for x in TE.finditer(labels)})
            category = "TE-related" if len(hits) >= 2 else "Possible_TE" if hits else "Review_non-TE"
            top = lambda key: "; ".join(f"{x} [{n}]" for x, n in counters[key].most_common(5)) or "NA"
            total = len(members[og])
            writer.writerow([og, total, annotated, f"{100*annotated/total:.1f}",
                             ";".join(hits) or "NA", category,
                             top("egg"), top("ipr"), top("pfam")])
    print(f"Selected {len(selected)} orthogroups; wrote {output}")


if __name__ == "__main__":
    main()
