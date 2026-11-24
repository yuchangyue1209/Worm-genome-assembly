#!/usr/bin/env python3
"""Plot a CAFE5 time tree with expanded/contracted/unchanged family counts.

Requires matplotlib; all input parsing uses the standard library. Internal
node IDs are read from a CAFE ASR tree and matched by descendant species,
never guessed from traversal order. No CAFE fitting or annotation is run.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from dataclasses import dataclass, field
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import re


ROOT = Path("/work/cyu/Schistocephalus_comparative_genomics")
RUN = "07_gene_family/cafe5/runs/screening8_16.5Ma_k3_v2/02_gamma_k3"
TIP_ORDER = ["E_granulosus", "E_multilocularis", "T_multiceps", "T_solium",
             "H_microstoma", "S_cotti", "S_solidus", "Spirometra_Aus1"]
LABELS = dict(zip(TIP_ORDER, ["E. granulosus", "E. multilocularis", "T. multiceps", "T. solium",
                            "H. microstoma", "S. cotti", "S. solidus", "Spirometra sp. Aus1"]))
# Okabe-Ito blue/vermillion, with a neutral background category.
# Keep +/- and named legend entries: meaning must not depend on color alone.
COLORS = {"Expanded": "#0072B2", "Contracted": "#D55E00", "Unchanged": "#D9D9D9"}
PALETTE_SOURCE = "https://jfly.uni-koeln.de/color/"


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float | None = None
    children: list = field(default_factory=list)
    parent: Node | None = None
    tips: frozenset = frozenset()
    node_id: str = ""
    depth: float = 0.0
    y: float = 0.0


def walk(node):
    yield node
    for child in node.children:
        yield from walk(child)


def parse_newick(text):
    text = re.sub(r"\[[^\]]*\]", "", text).strip()
    pos = 0
    def space():
        nonlocal pos
        while pos < len(text) and text[pos].isspace():
            pos += 1
    def token():
        nonlocal pos
        space()
        if pos < len(text) and text[pos] in "\"'":
            quote = text[pos]
            pos += 1
            chars = []
            while pos < len(text):
                if text[pos] == quote:
                    pos += 1
                    if pos < len(text) and text[pos] == quote:
                        chars.append(quote)
                        pos += 1
                        continue
                    return "".join(chars)
                chars.append(text[pos])
                pos += 1
            raise ValueError("Unclosed quoted Newick label")
        start = pos
        while pos < len(text) and text[pos] not in ",():;":
            pos += 1
        return text[start:pos].strip()
    def subtree():
        nonlocal pos
        space()
        node = Node()
        if pos < len(text) and text[pos] == "(":
            pos += 1
            node.children.append(subtree())
            space()
            while pos < len(text) and text[pos] == ",":
                pos += 1
                node.children.append(subtree())
                space()
            if pos >= len(text) or text[pos] != ")":
                raise ValueError("Newick missing closing parenthesis")
            pos += 1
        node.name = token()
        space()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            node.length = float(token())
        for child in node.children:
            child.parent = node
        return node
    root = subtree()
    space()
    if text[pos:] != ";":
        raise ValueError("Expected exactly one semicolon-terminated Newick tree")
    return root


def first_asr_tree(path):
    collecting, parts = False, []
    with path.open() as handle:
        for line in handle:
            if not collecting:
                match = re.match(r"\s*(?:tree|utree)\s+[^=]+?=\s*(.*)", line, re.I)
                if match:
                    collecting, line = True, match.group(1)
                elif line.lstrip().startswith("("):
                    collecting = True
                else:
                    continue
            parts.append(line.strip())
            if ";" in line:
                return parse_newick(" ".join(parts).split(";", 1)[0]+";")
    raise ValueError("No complete labelled CAFE tree found in "+str(path))


def clades(node, asr=False):
    if asr:
        match = re.fullmatch(r"(.*?)<(\d+)>\*?(?:_-?\d+)?", node.name)
        if not match:
            raise ValueError("Unrecognized ASR node label: "+node.name)
        prefix, node.node_id = match.groups()
        if not node.children:
            if not prefix:
                raise ValueError("ASR tip lacks a species name")
            node.tips = frozenset([prefix])
    elif not node.children:
        if not node.name:
            raise ValueError("Unnamed time-tree tip")
        node.tips = frozenset([node.name])
    if node.children:
        if len(node.children) != 2:
            raise ValueError("Expected a binary rooted tree")
        children = [clades(child, asr) for child in node.children]
        if children[0] & children[1]:
            raise ValueError("Duplicate tree species")
        node.tips = frozenset().union(*children)
    return node.tips


def map_nodes(dated, asr):
    clades(dated)
    clades(asr, True)
    lookup = {n.tips: n for n in walk(asr)}
    if len(lookup) != len(list(walk(asr))) or {n.tips for n in walk(dated)} != set(lookup):
        raise ValueError("ASR and supplied time tree differ in species, rooted topology or clades")
    ids = [n.node_id for n in walk(asr)]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate CAFE ASR node IDs")
    for node in walk(dated):
        source = lookup[node.tips]
        node.node_id = source.node_id
        if node.parent is not None:
            if node.length is None or not math.isfinite(node.length) or node.length <= 0:
                raise ValueError("Time-tree edge lengths must be finite and positive")
            if source.length is None or not math.isclose(node.length, source.length, rel_tol=3e-4, abs_tol=1e-4):
                raise ValueError("ASR branch length disagrees with supplied time tree at <"+node.node_id+">")
            node.depth = node.parent.depth+node.length
    depths = [n.depth for n in walk(dated) if not n.children]
    height = sum(depths)/len(depths)
    if height <= 0 or max(depths)-min(depths) > 1e-5*max(1.0, height):
        raise ValueError("Supplied time tree is not ultrametric")
    return height


def read_table(path):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [s.strip().lstrip("#") for s in next(reader, [])]
        while header and not header[-1]:
            header.pop()
        if not header or len(header) != len(set(header)):
            raise ValueError("Missing/duplicate header: "+str(path))
        rows = {}
        for i, cells in enumerate(reader, 2):
            if not cells or not any(s.strip() for s in cells):
                continue
            while len(cells) > len(header) and not cells[-1].strip():
                cells.pop()
            if len(cells) != len(header) or not cells[0].strip() or cells[0].strip() in rows:
                raise ValueError(f"Malformed/duplicate row {i}: {path}")
            rows[cells[0].strip()] = dict(zip(header, [s.strip() for s in cells]))
    return header, rows


def resolve(directory, stem):
    for ext in ("tab", "txt", "tsv"):
        path = directory/(stem+"."+ext)
        if path.is_file() and path.stat().st_size:
            return path
    raise FileNotFoundError("Missing nonempty "+str(directory/stem)+".tab/.txt/.tsv")


def column_ids(header):
    result = {}
    for column in header[1:]:
        match = re.fullmatch(r".*?<(\d+)>", column)
        if not match or match[1] in result:
            raise ValueError("Missing/duplicate explicit CAFE node ID: "+column)
        result[match[1]] = column
    return result


def calculate_stats(root, change_path, count_path, family_path=None, threshold=None):
    ch, changes = read_table(change_path)
    co, counts = read_table(count_path)
    if not changes or set(changes) != set(counts):
        raise ValueError("Count/change tables must contain the same nonempty family set")
    change_columns, count_columns = column_ids(ch), column_ids(co)
    all_ids = {n.node_id for n in walk(root)}
    branch_ids = all_ids-{root.node_id}
    if not branch_ids <= change_columns.keys() or not branch_ids <= count_columns.keys():
        raise ValueError("Reconstruction tables lack a mapped non-root branch")
    if set(change_columns)-all_ids or set(count_columns)-all_ids:
        raise ValueError("Reconstruction table contains nodes absent from the ASR tree")
    for node in walk(root):
        if node.children:
            continue
        expected = next(iter(node.tips))
        for col in (change_columns[node.node_id], count_columns[node.node_id]):
            if col != f"{expected}<{node.node_id}>":
                raise ValueError("Tip species/node ID mismatch: "+col)
    for row in counts.values():
        if any(int(row[col]) < 0 for col in count_columns.values()):
            raise ValueError("Negative reconstructed copy count")
    selected = set(changes)
    if threshold is not None:
        fh, families = read_table(family_path)
        pcol = next((s for s in fh if s.lower() == "pvalue"), None)
        if pcol is None or set(families) != set(changes):
            raise ValueError("Family p-values must cover the same reconstructed families")
        selected = set()
        for family, row in families.items():
            raw = row[pcol]
            if raw.upper() in {"NA", "NAN", "-NAN", "N/A", ""}:
                continue
            value = float(raw)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("Invalid family p-value: "+raw)
            if value < threshold:
                selected.add(family)
    if not selected:
        raise ValueError("No families in the selected display subset")
    stats = {}
    for node in walk(root):
        if node.parent is None:
            continue
        col, count_col = change_columns[node.node_id], count_columns[node.node_id]
        deltas = {f: int(row[col]) for f, row in changes.items()}
        parent_col = count_columns.get(node.parent.node_id)
        if parent_col:
            inconsistent = [f for f in changes if int(counts[f][count_col])-deltas[f] != int(counts[f][parent_col])]
            if inconsistent:
                raise ValueError("Count/change/parent reconstruction mismatch at "+col+": "+inconsistent[0])
        def summarize(ids):
            return dict(Expanded=sum(deltas[f] > 0 for f in ids), Contracted=sum(deltas[f] < 0 for f in ids),
                        Unchanged=sum(deltas[f] == 0 for f in ids))
        stats[node.node_id] = dict(Node_ID=node.node_id, Branch=col,
            Node_type="Internal" if node.children else "Tip", Descendant_species=",".join(sorted(node.tips)),
            Families=len(selected), **summarize(selected),
            All_expanded=summarize(changes)["Expanded"], All_contracted=summarize(changes)["Contracted"])
    return stats, len(changes), len(selected)


def check_clade_summary(path, stats, root_id=None):
    header, rows = read_table(path)
    if not {"Increase", "Decrease"} <= set(header):
        raise ValueError("Unexpected clade-results header")
    checked = set()
    for name, row in rows.items():
        match = re.fullmatch(r".*?<(\d+)>", name)
        if not match:
            raise ValueError("Unrecognized clade-results node: "+name)
        nid = match[1]
        if nid not in stats:
            if nid == root_id:
                continue  # Root has no incoming branch and is never plotted.
            raise ValueError("Unknown clade-summary node ID: "+nid)
        if (int(row["Increase"]), int(row["Decrease"])) != (stats[nid]["All_expanded"], stats[nid]["All_contracted"]):
            raise ValueError("Clade summary disagrees with family change signs at "+name)
        if nid in checked:
            raise ValueError("Duplicate clade-summary node ID")
        checked.add(nid)
    if checked != set(stats):
        raise ValueError("Clade summary lacks mapped branches")


def layout(root, order):
    rank = {sp: i for i, sp in enumerate(order)}
    for node in walk(root):
        node.children.sort(key=lambda child: min(rank[sp] for sp in child.tips))
    leaves = [n for n in walk(root) if not n.children]
    for i, node in enumerate(leaves):
        node.y = len(leaves)-1-i
    def position(node):
        if node.children:
            for child in node.children:
                position(child)
            node.y = sum(c.y for c in node.children)/len(node.children)
    position(root)
    return leaves


def draw(root, stats, height, order, args):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import DrawingArea, TextArea, HPacker, VPacker, AnnotationBbox
        from matplotlib.patches import Wedge, Patch
        from matplotlib.ticker import MaxNLocator
    except ImportError as exc:
        raise ValueError("matplotlib is required; use /home/cyu/.conda/envs/deeptmhmm2/bin/python") from exc
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "axes.unicode_minus": False})
    scale = args.font_scale
    leaves = layout(root, order)
    fig = plt.figure(figsize=(args.width, args.height), facecolor="white")
    ax = fig.add_axes([.065, .14, .60, .73])
    labels = fig.add_axes([.695, .14, .28, .73], sharey=ax)
    ax.set_ylim(-.85, len(leaves)+.15)
    ax.set_xlim(-height*.025, height*1.01)
    labels.set_xlim(0, 1)
    labels.axis("off")
    ax.set_yticks([])
    for name in ("left", "top", "right"):
        ax.spines[name].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.25)
    ages = [x for x in MaxNLocator(nbins=6).tick_values(0, height) if 0 <= x <= height]
    ax.set_xticks([height-x for x in ages], [f"{x:g}" for x in ages])
    ax.tick_params(axis="x", labelsize=14*scale, length=6, width=1.1, pad=7)
    ax.set_xlabel("Million years ago", fontsize=15*scale, labelpad=12)
    for node in walk(root):
        if node.children:
            ys = [c.y for c in node.children]
            ax.plot([node.depth]*2, [min(ys), max(ys)], color="#222222", lw=1.5, zorder=1)
        if node.parent is not None:
            ax.plot([node.parent.depth, node.depth], [node.y]*2, color="#222222", lw=1.5, zorder=1)
    def pie(values, diameter=32):
        diameter *= scale
        area = DrawingArea(diameter, diameter)
        radius = diameter/2
        total, theta = sum(values[x] for x in COLORS), 90.0
        for category in ("Expanded", "Contracted", "Unchanged"):
            angle = 360*values[category]/total
            if angle:
                area.add_artist(Wedge((radius, radius), radius, theta, theta+angle,
                    facecolor=COLORS[category], edgecolor="white", linewidth=.65))
            theta += angle
        return area
    def numbers(row, fontsize):
        return HPacker(children=[
            TextArea(f"+{row['Expanded']:,}", textprops=dict(color=COLORS["Expanded"], fontsize=fontsize)),
            TextArea(" / ", textprops=dict(color="#555555", fontsize=fontsize)),
            TextArea(f"-{row['Contracted']:,}", textprops=dict(color=COLORS["Contracted"], fontsize=fontsize))],
            align="center", pad=0, sep=0)
    for node in leaves:
        sp = next(iter(node.tips))
        focal = sp in {"S_cotti", "S_solidus"}
        color = "#151515"  # Focal species are bold, not change-category colors.
        if sp == "Spirometra_Aus1":
            name_box = HPacker(children=[
                TextArea("Spirometra", textprops=dict(fontsize=18*scale, fontstyle="italic", color=color)),
                TextArea(" sp. Aus1", textprops=dict(fontsize=18*scale, color=color))],
                align="center", pad=0, sep=0)
            labels.add_artist(AnnotationBbox(name_box, (0, node.y+.14), xycoords=labels.transData,
                box_alignment=(0, .5), frameon=False, pad=0))
        else:
            labels.text(0, node.y+.14, LABELS.get(sp, sp.replace("_", " ")),
                        fontsize=18*scale, fontstyle="italic", fontweight="bold" if focal else "normal",
                        color=color, va="center")
        labels.add_artist(AnnotationBbox(numbers(stats[node.node_id], 15*scale), (0, node.y-.25),
            xycoords=labels.transData, box_alignment=(0, .5), frameon=False, pad=0))
        labels.add_artist(AnnotationBbox(pie(stats[node.node_id], 34), (.94, node.y),
            xycoords=labels.transData, frameon=False, pad=0))
    handles = [Patch(facecolor=COLORS[x], label=x) for x in ("Expanded", "Contracted", "Unchanged")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.065, .965), frameon=False,
               ncol=3, fontsize=14*scale, handlelength=1.1, columnspacing=1.5)
    families = next(iter(stats.values()))["Families"]
    denominator_label = f"{families:,} families" if args.family_p_threshold is None else f"Family-wide p < {args.family_p_threshold:g}: {families:,} families"
    fig.text(.975, .943, denominator_label, ha="right", va="center", fontsize=12*scale, color="#555555")
    # Keep small internal branches at their real lengths. Move annotation boxes
    # rather than stretching branches; thin leaders preserve branch association.
    fig.canvas.draw()
    occupied, placements = [], []
    for node in sorted((n for n in walk(root) if n.children and n.parent is not None), key=lambda n: -n.y):
        row = stats[node.node_id]
        group = VPacker(children=[pie(row), numbers(row, 14*scale)], align="center", pad=0, sep=3)
        anchor = ((node.parent.depth+node.depth)/2, node.y)
        label = AnnotationBbox(group, anchor, xybox=(0, 5), xycoords=ax.transData,
            boxcoords="offset points", box_alignment=(.5, 0), frameon=False, pad=0,
            arrowprops=dict(arrowstyle="-", lw=.65, color="#888888"))
        ax.add_artist(label)
        offsets = [(0, 5), (-24, 5), (24, 5), (-48, 5), (48, 5), (0, -60),
                   (-24, -60), (24, -60), (-72, 5), (72, 5), (-48, -60), (48, -60)]
        for dx, dy in offsets:
            label.xybox = (dx*scale, dy*scale)
            fig.canvas.draw()
            bbox = group.get_window_extent(fig.canvas.get_renderer()).expanded(1.06, 1.07)
            bounds = ax.get_window_extent()
            if bounds.contains(bbox.x0, bbox.y0) and bounds.contains(bbox.x1, bbox.y1) and not any(bbox.overlaps(b) for b in occupied):
                occupied.append(bbox)
                placements.append(dict(Node_ID=node.node_id, offset_x_points=dx*scale, offset_y_points=dy*scale))
                break
        else:
            plt.close(fig)
            raise ValueError("Internal annotations overlap; increase --width/--height or reduce --font-scale")
    fig.canvas.draw()
    return fig, placements


def write_table(path, columns, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--tree", type=Path)
    parser.add_argument("--asr", type=Path, help="CAFE labelled ASR .tre/.nwk; used only for node mapping")
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--model", choices=["Gamma", "Base"], default="Gamma")
    parser.add_argument("--expected-species", type=int, default=8)
    parser.add_argument("--tip-order", nargs="+")
    parser.add_argument("--family-p-threshold", type=float, help="Optional exploratory family-wide subset, not branch significance")
    parser.add_argument("--font-scale", type=float, default=1.1)
    parser.add_argument("--width", type=float, default=17.5)
    parser.add_argument("--height", type=float, default=10.5)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Output exists; choose a new version: "+str(args.out))
    if any(not math.isfinite(x) or x <= 0 for x in (args.font_scale, args.width, args.height)) or args.dpi < 72:
        raise ValueError("Invalid font scale, figure dimensions or DPI")
    if args.family_p_threshold is not None and not 0 < args.family_p_threshold < 1:
        raise ValueError("Family p threshold must be between zero and one")
    results = args.results or args.root/RUN
    tree_path = args.tree or args.root/"07_gene_family/cafe5/time_trees/screening8_dated_16.5Ma.nwk"
    asr_path = args.asr or results/(args.model+"_asr.tre")
    if not asr_path.is_file():
        raise FileNotFoundError("CAFE node-labelled ASR tree required; supply --asr PATH. Do not infer internal IDs by order: "+str(asr_path))
    root, asr = parse_newick(tree_path.read_text()), first_asr_tree(asr_path)
    height = map_nodes(root, asr)
    if len(root.tips) != args.expected_species:
        raise ValueError("Unexpected number of time-tree species")
    order = args.tip_order or [x for x in TIP_ORDER if x in root.tips]+sorted(root.tips-set(TIP_ORDER))
    if len(order) != len(set(order)) or set(order) != root.tips:
        raise ValueError("Tip order must contain each tree species once")
    paths = dict(tree=tree_path, labelled_asr=asr_path,
                 change=resolve(results, args.model+"_change"), count=resolve(results, args.model+"_count"))
    if args.family_p_threshold is not None:
        paths["family_tests"] = resolve(results, args.model+"_family_results")
    stats, all_families, selected = calculate_stats(root, paths["change"], paths["count"],
        paths.get("family_tests"), args.family_p_threshold)
    try:
        clade_path = resolve(results, args.model+"_clade_results")
    except FileNotFoundError:
        clade_path = None
    if clade_path:
        check_clade_summary(clade_path, stats, root.node_id)
        paths["clade_summary"] = clade_path
    for node in walk(root):
        if node.node_id in stats:
            stats[node.node_id]["Node_age_Ma"] = max(0, height-node.depth)
    figure, placements = draw(root, stats, height, order, args)
    inputs = {key: dict(path=str(path.resolve()), sha256=sha256(path)) for key, path in paths.items()}
    args.out.mkdir(parents=True)
    for ext in ("pdf", "svg", "png"):
        figure.savefig(args.out/("cafe5_tree_pies."+ext), dpi=args.dpi, facecolor="white", bbox_inches="tight", pad_inches=.15)
    write_table(args.out/"branch_family_counts.tsv", ["Node_ID", "Branch", "Node_type", "Descendant_species",
        "Node_age_Ma", "Families", "Expanded", "Contracted", "Unchanged", "All_expanded", "All_contracted"], stats.values())
    write_table(args.out/"node_mapping.tsv", ["Node_ID", "Node_type", "Descendant_species", "Node_age_Ma", "Incoming_branch"],
        [dict(Node_ID=n.node_id, Node_type="Root" if n.parent is None else "Internal" if n.children else "Tip",
              Descendant_species=",".join(sorted(n.tips)), Node_age_Ma=max(0, height-n.depth),
              Incoming_branch="none_root" if n.parent is None else stats[n.node_id]["Branch"]) for n in walk(root)])
    (args.out/"figure_provenance.json").write_text(json.dumps(dict(
        created_UTC=dt.datetime.now(dt.timezone.utc).isoformat(), inputs=inputs, model=args.model,
        all_reconstructed_families=all_families, displayed_families=selected,
        family_p_threshold=args.family_p_threshold, source_files_modified=False,
        node_mapping="CAFE ASR explicit IDs matched by rooted descendant-species sets",
        clade_summary_crosscheck="passed" if clade_path else "not_available",
        branch_probabilities_used_as_pvalues=False, root_height_Ma=height,
        tip_order=[next(iter(n.tips)) for n in walk(root) if not n.children],
        font_scale=args.font_scale, figure_inches=[args.width, args.height], internal_annotation_offsets=placements,
        palette=COLORS, palette_source=PALETTE_SOURCE,
        redundant_encoding="named legend and signed +/- counts; focal species in bold black",
        interpretation="counts of families with reconstructed changes, not copy totals or branch significance",
        time_axis="supplied screening time tree; no independent age validation or confidence intervals"), indent=2)+"\n")
    (args.out/"interpretation_notes.txt").write_text(
        "Blue/vermillion/light-gray pies: expanded/contracted/unchanged families on the incoming branch.\n"
        "Red-green color-vision-friendly palette: Okabe-Ito blue #0072B2 and vermillion #D55E00, with neutral #D9D9D9. Signed +/- counts and named legend entries provide redundant encoding. Focal species are bold black.\n"
        "All pies use the displayed family subset as denominator. An absent family with no copy change is included as unchanged.\n"
        "Numbers count families with positive/negative reconstructed changes, not total gene copies gained/lost.\n"
        "Default: all reconstructed families, without significance filtering. No branch posterior is used as a p-value.\n"
        "Optional family-p subset is nominal family-wide screening plus branch change, not branch-specific significance or FDR. Reported p=0 is not an exact zero probability.\n"
        "Node mapping uses the first ASR tree only to recover stable explicit node IDs/topology; its single-family counts are never used as genome-wide totals.\n"
        "The root has no incoming branch and receives no expansion/contraction pie. No WGT, root gene count, age confidence interval or pooled-change pie is invented.\n"
        "Branch lengths come from the supplied time tree and are crosschecked against ASR lengths; labels may be moved with leaders but branches are never stretched.\n"
        "The default dated tree uses an analyst-selected 16.5-Ma secondary calibration. Its source IQ-TREE was the without32 sensitivity tree, whereas family counts use the selected baseline. These provenance limits remain applicable.\n")
    print("PASS:", args.out/"cafe5_tree_pies.pdf")
    print(f"Species={len(root.tips)}; reconstructed families={all_families}; displayed families={selected}")
    for row in stats.values():
        if row["Node_type"] == "Tip":
            print(row["Descendant_species"], f"+{row['Expanded']}/-{row['Contracted']}", sep="\t")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        raise SystemExit("STOPPED: "+str(exc))
