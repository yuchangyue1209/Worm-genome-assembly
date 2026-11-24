#!/usr/bin/env python3
"""Plot a clean large-font phylogeny and OrthoFinder summary for 11 cestodes.

The phylogeny is read from an IQ-TREE Newick file. Bar values are read from
the selected OrthoFinder result directory. The figure deliberately has no
overall title, legend, explanatory subtitle, footnote, or generated caption.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.ticker import MaxNLocator, PercentFormatter


LABELS = {
    "E_granulosus": "E. granulosus", "E_multilocularis": "E. multilocularis",
    "T_multiceps": "T. multiceps", "T_solium": "T. solium",
    "H_microstoma": "H. microstoma", "S_cotti": "S. cotti",
    "S_solidus": "S. solidus", "D_latus": "D. latus",
    "L_intestinalis": "L. intestinalis", "S_proliferum": "Sparganum proliferum",
    "Spirometra_Aus1": "Spirometra sp. Aus1",
}


def label_parts(species):
    """Italic biological name; roman rank abbreviation and isolate identifier."""
    if species == "Spirometra_Aus1":
        return [("Spirometra", "italic"), (" sp. Aus1", "normal")]
    return [(LABELS[species], "italic")]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    other = all_species - side
    a = tuple(sorted(side))
    b = tuple(sorted(other))
    return min(a, b)


def short_support(label: str) -> str:
    """Keep SH-aLRT/UFBoot and omit any appended concordance value."""
    fields = label.split("/")
    return "/".join(fields[:2]) if len(fields) >= 2 else label


def reroot_on_order_edge(original: Node, order_species: set[str]):
    nodes = list(walk(original))
    all_species = tips(original)
    adjacency = {node: [] for node in nodes}
    support = {}

    for node in nodes:
        for child in node.children:
            adjacency[node].append((child, child.length))
            adjacency[child].append((node, child.length))
            if child.children and child.name:
                support[split_key(tips(child), all_species)] = short_support(child.name)

    def component(start: Node, blocked: Node) -> set[str]:
        found = set()
        stack = [(start, blocked)]
        while stack:
            node, previous = stack.pop()
            if not node.children and len(adjacency[node]) == 1:
                found.add(node.name)
            for neighbour, _ in adjacency[node]:
                if neighbour is not previous:
                    stack.append((neighbour, node))
        return found

    cut = None
    seen = set()
    for left in nodes:
        for right, length in adjacency[left]:
            edge_id = frozenset((id(left), id(right)))
            if edge_id in seen:
                continue
            seen.add(edge_id)
            side = component(left, right)
            if side == order_species or all_species - side == order_species:
                cut = (left, right, length) if side == order_species else (right, left, length)
                break
        if cut:
            break
    if cut is None:
        raise ValueError("Could not find the edge separating Cyclophyllidea")

    cyclo_side, other_side, edge_length = cut

    def orient(node: Node, previous: Node, length: float) -> Node:
        return Node(
            name=node.name,
            length=length,
            children=[
                orient(neighbour, node, neighbour_length)
                for neighbour, neighbour_length in adjacency[node]
                if neighbour is not previous
            ],
        )

    root = Node(
        children=[
            orient(cyclo_side, other_side, edge_length / 2),
            orient(other_side, cyclo_side, edge_length / 2),
        ]
    )
    return root, support, all_species


def read_orthofinder_statistics(results: Path, expected_species: set[str]):
    stats_file = results / "Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv"
    duplication_file = results / "Comparative_Genomics_Statistics/Duplications_per_Species_Tree_Node.tsv"

    with stats_file.open() as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    species = rows[0][1:]
    if set(species) != expected_species:
        raise ValueError("IQ-TREE and OrthoFinder species sets differ")

    metrics = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            break
        if len(row) != len(species) + 1:
            raise ValueError(f"Malformed Statistics_PerSpecies row: {row[:2]}")
        metrics[row[0]] = {sp: float(value) for sp, value in zip(species, row[1:])}

    required = {
        "Number of genes",
        "Number of genes in orthogroups",
        "Percentage of genes in orthogroups",
        "Number of species-specific orthogroups",
    }
    missing = required - metrics.keys()
    if missing:
        raise ValueError(f"Missing OrthoFinder metrics: {sorted(missing)}")

    terminal_duplications = {}
    with duplication_file.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            node = row["Species Tree Node"]
            if node in expected_species:
                terminal_duplications[node] = int(row["Duplications (50% support)"])
    if set(terminal_duplications) != expected_species:
        raise ValueError("Terminal duplication rows are incomplete")
    return metrics, terminal_duplications


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True, help="IQ-TREE .treefile or concordance Newick")
    parser.add_argument("--results", type=Path, required=True, help="OrthoFinder Results_* directory")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--font-scale", type=float, default=1.50)
    parser.add_argument("--width", type=float, default=18.5, help="Figure width in inches")
    parser.add_argument("--height", type=float, default=9.6, help="Figure height in inches")
    parser.add_argument("--no-branch-breaks", action="store_true", help="Draw all branches at their full lengths")
    parser.add_argument("--branch-break-threshold", type=float, default=.20)
    parser.add_argument("--compressed-length", type=float, default=.085)
    args = parser.parse_args()
    if any(not math.isfinite(value) or value <= 0 for value in
           (args.font_scale, args.width, args.height, args.branch_break_threshold, args.compressed_length)):
        raise ValueError("Font scale and figure dimensions must be finite and positive")
    if args.compressed_length >= args.branch_break_threshold:
        raise ValueError("Compressed display length must be smaller than the break threshold")
    if args.out.exists():
        raise FileExistsError(f"Output exists; choose a new version: {args.out}")
    fs = lambda value: value * args.font_scale

    original = parse_newick(args.tree.read_text())
    cyclophyllidea = {
        "E_granulosus", "E_multilocularis", "T_multiceps", "T_solium", "H_microstoma"
    }
    root, support, all_species = reroot_on_order_edge(original, cyclophyllidea)
    if len(all_species) != 11:
        raise ValueError(f"Expected 11 species, found {len(all_species)}")

    preference = [
        "E_granulosus", "E_multilocularis", "T_multiceps", "T_solium", "H_microstoma",
        "S_cotti", "S_solidus", "D_latus", "L_intestinalis",
        "S_proliferum", "Spirometra_Aus1",
    ]
    rank = {species: index for index, species in enumerate(preference)}

    def rotate(node: Node):
        for child in node.children:
            rotate(child)
        node.children.sort(key=lambda child: min(rank[species] for species in tips(child)))

    rotate(root)
    metrics, duplications = read_orthofinder_statistics(args.results, all_species)

    colours = {
        species: ("#008C7A" if species in cyclophyllidea else "#8B4F83")
        for species in all_species
    }
    colours["S_cotti"] = "#DF8500"
    colours["S_solidus"] = "#C9443B"

    leaves = [node for node in walk(root) if not node.children]
    order = [node.name for node in leaves]
    y = {species: index for index, species in enumerate(order)}
    xpos, ypos, break_rows = {}, {}, []

    def layout(node: Node, x: float = 0.0):
        xpos[node] = x
        for child in node.children:
            shown = child.length
            if not args.no_branch_breaks and child.length > args.branch_break_threshold:
                shown = args.compressed_length
                break_rows.append(dict(descendants=",".join(sorted(tips(child))),
                                       actual_length=child.length, displayed_length=shown))
            layout(child, x + shown)
        ypos[node] = (
            sum(ypos[child] for child in node.children) / len(node.children)
            if node.children else y[node.name]
        )

    layout(root)
    max_tree_x = max(xpos.values())
    guide_end = max_tree_x * 1.02
    label_x = max_tree_x * 1.04

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": fs(11),
        "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8,
    })
    fig = plt.figure(figsize=(args.width, args.height), facecolor="white")
    grid = fig.add_gridspec(
        1, 4, width_ratios=[6.0, 2.8, 3.1, 3.0],
        left=0.028, right=0.988, bottom=0.12, top=0.86, wspace=0.22,
    )
    tree_ax = fig.add_subplot(grid[0])
    bar_axes = [fig.add_subplot(grid[index]) for index in (1, 2, 3)]
    all_axes = [tree_ax, *bar_axes]
    for axis in all_axes:
        # Leave a full row below the last species for the phylogenetic scale bar.
        axis.set_ylim(len(order) + 0.65, -0.65)
        for species in ("S_cotti", "S_solidus"):
            axis.axhspan(y[species] - 0.43, y[species] + 0.43,
                        color=colours[species], alpha=0.075, lw=0, zorder=0)

    name_artists = []
    for node in walk(root):
        if node.children:
            tree_ax.plot(
                [xpos[node], xpos[node]],
                [min(ypos[child] for child in node.children), max(ypos[child] for child in node.children)],
                color="#33383D", lw=1.8, zorder=2,
            )
            for child in node.children:
                tree_ax.plot(
                    [xpos[node], xpos[child]], [ypos[child], ypos[child]],
                    color="#33383D", lw=1.8, zorder=2,
                )
                if not args.no_branch_breaks and child.length > args.branch_break_threshold:
                    midpoint = (xpos[node] + xpos[child]) / 2
                    gap = args.compressed_length * .075
                    tree_ax.plot([midpoint-gap, midpoint+gap], [ypos[child]]*2,
                                 color="white", lw=6, solid_capstyle="butt", zorder=3)
                    for delta in [-gap*.55, gap*.55]:
                        tree_ax.plot([midpoint+delta-gap*.35, midpoint+delta+gap*.35],
                                     [ypos[child]+.12, ypos[child]-.12],
                                     color="#33383D", lw=1.5, zorder=4)
                if child.children:
                    value = support.get(split_key(tips(child), all_species), "")
                    if value:
                        tree_ax.text(
                            xpos[child] - max_tree_x * 0.008, ypos[child] - 0.18, value,
                            ha="right", va="bottom", fontsize=fs(10.5), color="#4A4F54",
                            bbox=dict(facecolor="white", edgecolor="none", pad=0.3), zorder=4,
                        )
        else:
            tree_ax.plot([xpos[node], guide_end], [ypos[node]]*2,
                         color="#C5C9CD", lw=.8, ls=(0, (2, 3)), zorder=1)
            parts = [TextArea(text, textprops=dict(fontsize=fs(14.0), fontstyle=style,
                     fontweight="bold" if node.name in {"S_cotti", "S_solidus"} else "normal",
                     color=colours[node.name])) for text, style in label_parts(node.name)]
            box = HPacker(children=parts, align="center", pad=0, sep=0)
            artist = AnnotationBbox(box, (label_x, ypos[node]),
                xycoords=tree_ax.transData, box_alignment=(0, .5), frameon=False, pad=0,
                annotation_clip=False, zorder=5)
            tree_ax.add_artist(artist)
            name_artists.append(artist)

    scale = 0.05
    scale_x = 0.01
    scale_y = len(order) + 0.10
    tree_ax.plot([scale_x, scale_x + scale], [scale_y, scale_y], color="#252525", lw=1.8)
    tree_ax.text(scale_x + scale / 2, scale_y + 0.27, "0.05 substitutions/site",
                 ha="center", va="top", fontsize=fs(10.5))
    tree_ax.set_xlim(-max_tree_x * 0.015, max_tree_x * 1.75)
    tree_ax.axis("off")
    tree_ax.text(0.10, 1.055, "Species phylogeny", transform=tree_ax.transAxes,
                 fontsize=fs(17), fontweight="bold")

    assigned_percent = metrics["Percentage of genes in orthogroups"]
    assigned_count = metrics["Number of genes in orthogroups"]
    species_specific = metrics["Number of species-specific orthogroups"]
    panel_values = [
        [duplications[species] for species in order],
        [assigned_percent[species] for species in order],
        [species_specific[species] for species in order],
    ]
    panel_titles = [
        "Terminal gene\nduplications", "Proteins in\northogroups",
        "Species-specific\northogroups",
    ]
    xlabels = ["Events", "Assigned proteins (%)", "Orthogroups"]

    for panel_index, axis in enumerate(bar_axes):
        values = panel_values[panel_index]
        xmax = (100.0 if panel_index == 1 else max(values))
        headroom = 1.17 if panel_index == 1 else 1.23
        axis.barh(
            range(len(order)), values, height=0.66,
            color=[colours[species] for species in order], edgecolor="none", zorder=3,
        )
        axis.set_xlim(0, xmax * headroom)
        axis.set_yticks([])
        axis.xaxis.grid(True, color="#E4E7EA", lw=0.8, zorder=1)
        axis.set_axisbelow(True)
        for side in ("top", "right", "left"):
            axis.spines[side].set_visible(False)
        axis.spines["bottom"].set_color("#AEB4B9")
        axis.tick_params(axis="x", labelsize=fs(10.5), length=3, color="#8E959B")
        axis.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=panel_index != 1))
        if panel_index == 1:
            axis.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        axis.set_xlabel(xlabels[panel_index], fontsize=fs(11.5), labelpad=8)
        axis.text(0, 1.055, panel_titles[panel_index], transform=axis.transAxes,
                  fontsize=fs(16), fontweight="bold", linespacing=1.02)

        for row_index, (species, value) in enumerate(zip(order, values)):
            label = (
                f"{int(assigned_count[species]):,}" if panel_index == 1
                else f"{int(value):,}"
            )
            axis.text(
                value + xmax * 0.025, row_index, label,
                ha="left", va="center", fontsize=fs(11.2),
                fontweight="bold" if species in {"S_cotti", "S_solidus"} else "normal",
                color="#30343A", zorder=5,
            )

    # Expand only blank label space. Shortened branches are explicitly marked //;
    # the 0.05 scale applies to all unbroken segments, never the dashed guides.
    for _ in range(12):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        overflow = max(artist.get_window_extent(renderer).x1 for artist in name_artists) - tree_ax.bbox.x1
        if overflow <= 0:
            break
        left, right = tree_ax.get_xlim()
        tree_ax.set_xlim(left, right + (right-left)*(overflow+12)/tree_ax.bbox.width)
    else:
        raise ValueError("Species labels do not fit; increase --width or reduce --font-scale")

    args.out.mkdir(parents=True, exist_ok=False)
    stem = args.out / "orthology_summary_clean"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=300, facecolor="white",
                    bbox_inches="tight", pad_inches=0.08)
    fig.savefig(args.out / "orthology_summary_clean_preview.png", dpi=120,
                facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    source_data = []
    for species in order:
        source_data.append(dict(species=species, display_name=LABELS[species],
            input_proteins=int(metrics["Number of genes"][species]),
            terminal_duplications_support_ge_0_5=duplications[species],
            assigned_proteins=int(assigned_count[species]), assigned_percent=assigned_percent[species],
            species_specific_orthogroups=int(species_specific[species])))
    with (args.out / "orthology_summary_clean.source_data.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, source_data[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(source_data)
    with (args.out / "orthology_summary_clean.branch_breaks.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, ["descendants", "actual_length", "displayed_length"], delimiter="\t")
        writer.writeheader()
        writer.writerows(break_rows)
    audit = {
        "tree": str(args.tree.resolve()),
        "orthofinder_results": str(args.results.resolve()),
        "species": order,
        "font_scale": args.font_scale,
        "figure_inches": [args.width, args.height],
        "label_parts": {species: label_parts(species) for species in order},
        "source_sha256": {str(path.resolve()): sha256(path) for path in [args.tree,
            args.results / "Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv",
            args.results / "Comparative_Genomics_Statistics/Duplications_per_Species_Tree_Node.tsv"]},
        "source_files_modified": False,
        "source_tree_branch_lengths_modified": False,
        "displayed_branches_compressed": bool(break_rows),
        "branch_break_threshold": args.branch_break_threshold,
        "compressed_length": args.compressed_length,
        "branch_breaks": break_rows,
        "scale_applies_to": "unbroken_segments_only" if break_rows else "all_branches",
        "tip_extensions": "dashed alignment guides, not evolutionary branch lengths",
        "outputs": [str(stem.with_suffix(f".{extension}")) for extension in ("png", "pdf", "svg")],
    }
    (args.out / "orthology_summary_clean.audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
