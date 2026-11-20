#!/usr/bin/env python3
"""Classify per-locus trees among the three competing diphyllobothriid roots."""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from Bio import Phylo


def monophyletic(tree, taxa: set[str]) -> bool:
    clade = tree.common_ancestor(list(taxa))
    observed = {tip.name for tip in clade.get_terminals()}
    return observed == taxa


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", type=Path, required=True,
                        help="One Newick gene tree per line")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    all_taxa = {"S_solidus", "S_cotti", "S_proliferum", "Spirometra_Aus1",
                "D_latus", "L_intestinalis"}
    candidates = {
        "Sparganum_Spirometra_basal": {"S_solidus", "S_cotti", "D_latus", "L_intestinalis"},
        "Schistocephalus_basal": {"S_proliferum", "Spirometra_Aus1", "D_latus", "L_intestinalis"},
        "D_latus_Ligula_basal": {"S_solidus", "S_cotti", "S_proliferum", "Spirometra_Aus1"},
    }
    counts = Counter()
    total = 0
    for tree in Phylo.parse(str(args.trees), "newick"):
        names = {tip.name for tip in tree.get_terminals()}
        if not all_taxa <= names:
            counts["Other_or_missing_taxa"] += 1
            total += 1
            continue
        supported = [label for label, ingroup in candidates.items() if monophyletic(tree, ingroup)]
        counts[supported[0] if len(supported) == 1 else "Other_or_unresolved"] += 1
        total += 1
    if total == 0:
        raise ValueError("No trees were read")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    order = list(candidates) + ["Other_or_unresolved", "Other_or_missing_taxa"]
    with args.output.open("w") as handle:
        handle.write("Topology\tGene_trees\tPercent\n")
        for label in order:
            if counts[label]:
                handle.write(f"{label}\t{counts[label]}\t{100 * counts[label] / total:.2f}\n")
    print(f"Wrote {args.output}; trees={total}")


if __name__ == "__main__":
    main()
