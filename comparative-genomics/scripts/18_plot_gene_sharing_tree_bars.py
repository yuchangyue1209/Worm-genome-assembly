#!/usr/bin/env python3
"""Plot an IQ-TREE phylogeny beside mutually exclusive gene-sharing classes.

Inputs are a matching IQ-TREE Newick and OrthoFinder Results_* directory.
Every representative gene/protein is counted exactly once. Unassigned genes
are combined with genes from species-specific orthogroups in the
"Species-specific / singleton" class.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)


def parse_newick(text: str) -> Node:
    text = text.strip()
    i = 0

    def token() -> str:
        nonlocal i
        start = i
        while i < len(text) and text[i] not in ",():;":
            i += 1
        return text[start:i].strip()

    def subtree() -> Node:
        nonlocal i
        node = Node()
        if i < len(text) and text[i] == "(":
            i += 1
            node.children.append(subtree())
            while text[i] == ",":
                i += 1
                node.children.append(subtree())
            if text[i] != ")":
                raise ValueError("Invalid Newick: missing ')'")
            i += 1
        node.name = token()
        if i < len(text) and text[i] == ":":
            i += 1
            node.length = float(token())
        return node

    root = subtree()
    if text[i:] != ";":
        raise ValueError(f"Unexpected Newick suffix: {text[i:]!r}")
    return root


def walk(node: Node):
    yield node
    for child in node.children:
        yield from walk(child)


def tips(node: Node) -> set[str]:
    if not node.children:
        return {node.name}
    return set().union(*(tips(child) for child in node.children))


def split_key(side: set[str], all_species: set[str]):
    a = tuple(sorted(side))
    b = tuple(sorted(all_species - side))
    return min(a, b)


def support_label(label: str) -> str:
    fields = label.split("/")
    return "/".join(fields[:2]) if len(fields) >= 2 else label


def reroot_between_orders(original: Node, cyclophyllidea: set[str]):
    nodes = list(walk(original))
    all_species = tips(original)
    adjacency = {node: [] for node in nodes}
    supports = {}

    for node in nodes:
        for child in node.children:
            adjacency[node].append((child, child.length))
            adjacency[child].append((node, child.length))
            if child.children and child.name:
                supports[split_key(tips(child), all_species)] = support_label(child.name)

    def component(start: Node, blocked: Node):
        result = set()
        stack = [(start, blocked)]
        while stack:
            node, previous = stack.pop()
            if not node.children and len(adjacency[node]) == 1:
                result.add(node.name)
            for neighbour, _ in adjacency[node]:
                if neighbour is not previous:
                    stack.append((neighbour, node))
        return result

    cut = None
    seen = set()
    for left in nodes:
        for right, length in adjacency[left]:
            edge = tuple(sorted((id(left), id(right))))
            if edge in seen:
                continue
            seen.add(edge)
            side = component(left, right)
            if side == cyclophyllidea:
                cut = left, right, length
                break
            if all_species - side == cyclophyllidea:
                cut = right, left, length
                break
        if cut:
            break
    if cut is None:
        raise ValueError("No edge separates the sampled Cyclophyllidea")

    cyclo_node, diphy_node, edge_length = cut

    def orient(node: Node, previous: Node, length: float):
        return Node(
            name=node.name,
            length=length,
            children=[
                orient(neighbour, node, neighbour_length)
                for neighbour, neighbour_length in adjacency[node]
                if neighbour is not previous
            ],
        )

    # Diphyllobothriidea is drawn above Cyclophyllidea, matching the reference style.
    root = Node(children=[
        orient(diphy_node, cyclo_node, edge_length / 2),
        orient(cyclo_node, diphy_node, edge_length / 2),
    ])
    return root, supports, all_species


def read_species_totals(results: Path, expected_species: set[str]):
    path = results / "Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv"
    with path.open() as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    species = rows[0][1:]
    if set(species) != expected_species:
        raise ValueError("Tree and OrthoFinder species sets differ")
    metrics = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            break
        if len(row) != len(species) + 1:
            raise ValueError(f"Malformed per-species row: {row[:2]}")
        metrics[row[0]] = {sp: int(float(value)) for sp, value in zip(species, row[1:])}
    for required in ("Number of genes", "Number of genes in orthogroups", "Number of unassigned genes"):
        if required not in metrics:
            raise ValueError(f"Missing metric: {required}")
    return species, metrics


def classify_gene_counts(results: Path, species: list[str], metrics: dict):
    categories = [
        "Shared by all (multi-copy)",
        "Shared by all (one-to-one)",
        "Genus-specific",
        "Species-specific / singleton",
        "Others",
    ]
    counts = {sp: Counter({category: 0 for category in categories}) for sp in species}
    assigned_from_gene_count = Counter()
    genera = {
        "Schistocephalus": {"S_cotti", "S_solidus"},
        "Echinococcus": {"E_granulosus", "E_multilocularis"},
        "Taenia": {"T_solium", "T_multiceps"},
    }

    table = results / "Orthogroups/Orthogroups.GeneCount.tsv"
    orthogroups_by_class = Counter()
    with table.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or []) - {"Orthogroup", "Total"} != set(species):
            raise ValueError("Unexpected Orthogroups.GeneCount.tsv header")
        for row in reader:
            values = {sp: int(row[sp]) for sp in species}
            if sum(values.values()) != int(row["Total"]):
                raise ValueError(f"Gene-count total mismatch in {row['Orthogroup']}")
            present = {sp for sp, value in values.items() if value > 0}

            if present == set(species):
                category = (
                    "Shared by all (one-to-one)"
                    if all(value == 1 for value in values.values())
                    else "Shared by all (multi-copy)"
                )
            elif len(present) == 1:
                category = "Species-specific / singleton"
            elif any(len(present) >= 2 and present <= members for members in genera.values()):
                category = "Genus-specific"
            else:
                category = "Others"

            orthogroups_by_class[category] += 1
            for sp, value in values.items():
                counts[sp][category] += value
                assigned_from_gene_count[sp] += value

    # Some OrthoFinder releases contain small disagreements between
    # Orthogroups.GeneCount.tsv and Statistics_PerSpecies.tsv. Use the actual
    # GeneCount memberships for classification, then calculate the residual
    # against the input total so every representative gene is counted once.
    reconciliation = {}
    for sp in species:
        total = metrics["Number of genes"][sp]
        residual = total - assigned_from_gene_count[sp]
        if residual < 0:
            raise ValueError(
                f"GeneCount memberships exceed input genes for {sp}: "
                f"{assigned_from_gene_count[sp]} > {total}"
            )
        counts[sp]["Species-specific / singleton"] += residual
        observed = sum(counts[sp].values())
        expected = total
        if observed != expected:
            raise ValueError(f"Category total mismatch for {sp}: {observed} != {expected}")
        reconciliation[sp] = {
            "input_genes": total,
            "assigned_from_gene_count": assigned_from_gene_count[sp],
            "assigned_from_statistics": metrics["Number of genes in orthogroups"][sp],
            "gene_count_minus_statistics": (
                assigned_from_gene_count[sp] - metrics["Number of genes in orthogroups"][sp]
            ),
            "singleton_residual_used": residual,
            "unassigned_from_statistics": metrics["Number of unassigned genes"][sp],
        }
    return categories, counts, orthogroups_by_class, reconciliation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--font-scale", type=float, default=1.15)
    args = parser.parse_args()
    if args.font_scale <= 0:
        raise ValueError("--font-scale must be greater than zero")
    args.out.mkdir(parents=True, exist_ok=True)
    fs = lambda value: value * args.font_scale

    cyclophyllidea = {
        "E_granulosus", "E_multilocularis", "H_microstoma", "T_multiceps", "T_solium"
    }
    original = parse_newick(args.tree.read_text())
    root, supports, all_species = reroot_between_orders(original, cyclophyllidea)
    if len(all_species) != 11:
        raise ValueError(f"Expected 11 species, found {len(all_species)}")

    preference = [
        "S_proliferum", "Spirometra_Aus1", "S_cotti", "S_solidus",
        "D_latus", "L_intestinalis", "H_microstoma", "T_multiceps",
        "T_solium", "E_granulosus", "E_multilocularis",
    ]
    rank = {species: index for index, species in enumerate(preference)}

    def rotate(node: Node):
        for child in node.children:
            rotate(child)
        node.children.sort(key=lambda child: min(rank[species] for species in tips(child)))

    rotate(root)
    species, metrics = read_species_totals(args.results, all_species)
    categories, class_counts, og_class_counts, reconciliation = classify_gene_counts(
        args.results, species, metrics
    )

    display = {
        "S_proliferum": "Sparganum proliferum", "Spirometra_Aus1": "Spirometra sp. Aus1",
        "S_cotti": "S. cotti", "S_solidus": "S. solidus",
        "D_latus": "D. latus", "L_intestinalis": "L. intestinalis",
        "H_microstoma": "H. microstoma", "T_multiceps": "T. multiceps",
        "T_solium": "T. solium", "E_granulosus": "E. granulosus",
        "E_multilocularis": "E. multilocularis",
    }
    order = [node.name for node in walk(root) if not node.children]
    y = {sp: index for index, sp in enumerate(order)}

    xpos, ypos = {}, {}

    def layout(node: Node, x: float = 0.0):
        xpos[node] = x
        for child in node.children:
            layout(child, x + child.length)
        ypos[node] = (
            sum(ypos[child] for child in node.children) / len(node.children)
            if node.children else y[node.name]
        )

    layout(root)
    max_tree_x = max(xpos.values())

    colours = {
        "Shared by all (multi-copy)": "#463078",
        "Shared by all (one-to-one)": "#83BFA2",
        "Genus-specific": "#FFF0A6",
        "Species-specific / singleton": "#F68F89",
        "Others": "#9B9B9B",
    }
    diphy_background = "#EAF8F7"
    cyclo_background = "#FFF8E8"

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": fs(11),
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(21.5, 10.5), facecolor="white")
    grid = fig.add_gridspec(
        1, 2, width_ratios=[1.42, 1.0], left=0.025, right=0.985,
        bottom=0.105, top=0.96, wspace=0.035,
    )
    tree_ax = fig.add_subplot(grid[0])
    bar_ax = fig.add_subplot(grid[1])
    for axis in (tree_ax, bar_ax):
        # Leave a full row below the last species for the phylogenetic scale bar.
        axis.set_ylim(len(order) + 0.65, -0.65)
        diphy_rows = [y[sp] for sp in order if sp not in cyclophyllidea]
        cyclo_rows = [y[sp] for sp in order if sp in cyclophyllidea]
        axis.axhspan(min(diphy_rows) - 0.55, max(diphy_rows) + 0.55,
                    color=diphy_background, lw=0, zorder=0)
        axis.axhspan(min(cyclo_rows) - 0.55, max(cyclo_rows) + 0.55,
                    color=cyclo_background, lw=0, zorder=0)

    for node in walk(root):
        if node.children:
            tree_ax.plot(
                [xpos[node], xpos[node]],
                [min(ypos[child] for child in node.children), max(ypos[child] for child in node.children)],
                color="black", lw=1.9, zorder=2,
            )
            for child in node.children:
                tree_ax.plot(
                    [xpos[node], xpos[child]], [ypos[child], ypos[child]],
                    color="black", lw=1.9, zorder=2,
                )
                if child.children:
                    label = supports.get(split_key(tips(child), all_species), "")
                    if label:
                        tree_ax.text(
                            xpos[child] + max_tree_x * 0.006, ypos[child] - 0.17, label,
                            ha="left", va="bottom", fontsize=fs(9.5), color="#202020",
                            bbox=dict(facecolor="white", alpha=0.78, edgecolor="none", pad=0.2),
                            zorder=4,
                        )
        else:
            colour = "#1F1F1F"
            if node.name == "S_cotti":
                colour = "#D97800"
            elif node.name == "S_solidus":
                colour = "#C23F38"
            tree_ax.plot(
                [xpos[node], max_tree_x * 1.20], [ypos[node], ypos[node]],
                color="#555555", lw=0.9, ls=(0, (3, 3)), zorder=1,
            )
            tree_ax.text(
                max_tree_x * 1.23, ypos[node], display[node.name],
                ha="left", va="center", fontsize=fs(13), fontstyle="italic",
                fontweight="bold" if node.name in {"S_cotti", "S_solidus"} else "normal",
                color=colour, zorder=5,
            )

    diphy_rows = [y[sp] for sp in order if sp not in cyclophyllidea]
    cyclo_rows = [y[sp] for sp in order if sp in cyclophyllidea]
    tree_ax.text(-max_tree_x * 0.12, sum(diphy_rows) / len(diphy_rows),
                 "Diphyllobothriidea", ha="right", va="center", fontsize=fs(13))
    tree_ax.text(-max_tree_x * 0.12, sum(cyclo_rows) / len(cyclo_rows),
                 "Cyclophyllidea", ha="right", va="center", fontsize=fs(13))

    scale = 0.05
    scale_x = max_tree_x * 0.05
    scale_y = len(order) + 0.10
    tree_ax.plot([scale_x, scale_x + scale], [scale_y, scale_y], color="black", lw=2)
    tree_ax.text(scale_x + scale / 2, scale_y + 0.27, "0.05",
                 ha="center", va="top", fontsize=fs(11))
    tree_ax.set_xlim(-max_tree_x * 0.38, max_tree_x * 2.40)
    tree_ax.axis("off")

    left = [0.0] * len(order)
    for category in categories:
        values = [class_counts[sp][category] for sp in order]
        bar_ax.barh(
            range(len(order)), values, left=left, height=0.62,
            color=colours[category], edgecolor="none", label=category, zorder=2,
        )
        left = [a + b for a, b in zip(left, values)]

    max_total = max(metrics["Number of genes"].values())
    bar_ax.set_xlim(0, max_total * 1.10)
    bar_ax.set_yticks([])
    bar_ax.set_xlabel("Number of genes", fontsize=fs(13), labelpad=10)
    bar_ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    bar_ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    bar_ax.tick_params(axis="x", labelsize=fs(10.5), direction="in", length=8, width=1.5)
    for side in ("top", "right", "left"):
        bar_ax.spines[side].set_visible(False)
    bar_ax.spines["bottom"].set_linewidth(1.6)

    legend_labels = [
        ("Others", "Others"),
        ("Species-specific / singleton", "Species-specific / singleton"),
        ("Genus-specific", "Genus-specific (≥2 sampled species)"),
        ("Shared by all (one-to-one)", "Shared by all (one-to-one)"),
        ("Shared by all (multi-copy)", "Shared by all (multi/one-to-multi)"),
    ]
    handles = [Patch(facecolor=colours[key], edgecolor="none", label=label)
               for key, label in legend_labels]
    bar_ax.legend(
        handles=handles, loc="center right", bbox_to_anchor=(0.99, 0.38),
        frameon=False, fontsize=fs(11.5), handlelength=1.4, handleheight=1.4,
        labelspacing=1.05,
    )

    stem = args.out / "gene_sharing_tree_bars"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=300, facecolor="white",
                    bbox_inches="tight", pad_inches=0.04)
    fig.savefig(args.out / "gene_sharing_tree_bars_preview.png", dpi=120,
                facecolor="white", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    with (args.out / "gene_sharing_counts.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Species", *categories, "Total"])
        for sp in order:
            values = [class_counts[sp][category] for category in categories]
            writer.writerow([sp, *values, sum(values)])

    audit = {
        "tree": str(args.tree.resolve()),
        "orthofinder_results": str(args.results.resolve()),
        "species": order,
        "orthogroups_by_class": dict(og_class_counts),
        "category_totals_equal_input_gene_counts": True,
        "gene_count_statistics_reconciliation": reconciliation,
        "genus_specific_definition": (
            "Orthogroup present in at least two sampled members of Schistocephalus, "
            "Echinococcus or Taenia and absent outside that genus"
        ),
    }
    (args.out / "gene_sharing_tree_bars.audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
