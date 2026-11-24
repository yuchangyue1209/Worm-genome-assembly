#!/usr/bin/env python3
"""Validate that a CAFE5 tree is binary, ultrametric and matches its table."""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Node:
    name: str = ""
    length: float | None = None
    children: list["Node"] = field(default_factory=list)


def parse_newick(text: str) -> Node:
    text = text.strip()
    position = 0

    def token() -> str:
        nonlocal position
        start = position
        while position < len(text) and text[position] not in ",():;":
            position += 1
        return text[start:position].strip().strip("'\"")

    def subtree() -> Node:
        nonlocal position
        node = Node()
        if position < len(text) and text[position] == "(":
            position += 1
            node.children.append(subtree())
            while position < len(text) and text[position] == ",":
                position += 1
                node.children.append(subtree())
            if position >= len(text) or text[position] != ")":
                raise ValueError("Invalid Newick: missing closing parenthesis")
            position += 1
        node.name = token()
        if position < len(text) and text[position] == ":":
            position += 1
            raw = token()
            if not raw:
                raise ValueError("Missing branch length")
            node.length = float(raw)
        return node

    root = subtree()
    if text[position:] != ";":
        raise ValueError(f"Unexpected Newick suffix: {text[position:]!r}")
    return root


def walk(node: Node):
    yield node
    for child in node.children:
        yield from walk(child)


def tip_distances(node: Node, distance: float = 0.0) -> dict[str, float]:
    if not node.children:
        if not node.name:
            raise ValueError("Unnamed tip")
        return {node.name: distance}
    result: dict[str, float] = {}
    for child in node.children:
        if child.length is None:
            raise ValueError(f"Missing branch length leading to {child.name or 'internal node'}")
        if not math.isfinite(child.length) or child.length < 0:
            raise ValueError("Branch lengths must be finite and non-negative")
        for name, value in tip_distances(child, distance + child.length).items():
            if name in result:
                raise ValueError(f"Duplicate tree tip: {name}")
            result[name] = value
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--relative-tolerance", type=float, default=1e-5)
    args = parser.parse_args()

    root = parse_newick(args.tree.read_text())
    nonbinary = [node for node in walk(root) if node.children and len(node.children) != 2]
    if nonbinary:
        raise ValueError(f"Tree is not fully binary: {len(nonbinary)} non-bifurcating nodes")

    distances = tip_distances(root)
    with args.families.open(newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) < 4 or header[:2] != ["Desc", "Family ID"]:
        raise ValueError("Unexpected CAFE5 family-table header")
    species = header[2:]
    if set(species) != set(distances):
        raise ValueError(
            f"Tree/table species mismatch; tree_only={sorted(set(distances)-set(species))}, "
            f"table_only={sorted(set(species)-set(distances))}"
        )

    minimum, maximum = min(distances.values()), max(distances.values())
    scale = max(maximum, 1.0)
    difference = maximum - minimum
    if difference > args.relative_tolerance * scale:
        raise ValueError(
            f"Tree is not ultrametric: root-to-tip range {minimum:.10g}-{maximum:.10g}"
        )

    print(f"PASS: binary ultrametric tree with {len(distances)} tips")
    print(f"Root-to-tip time: {sum(distances.values())/len(distances):.10g}")
    print("NOTE: numeric ultrametricity does not prove that branch lengths use real time units.")


if __name__ == "__main__":
    main()
