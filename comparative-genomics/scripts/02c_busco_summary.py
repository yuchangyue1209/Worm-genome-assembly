#!/usr/bin/env python3
"""Combine protein-mode BUSCO summaries for the configured species set."""

import argparse
import re
from pathlib import Path


PATTERN = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],"
    r"F:([\d.]+)%,M:([\d.]+)%,n:(\d+)"
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--busco-root", required=True, type=Path)
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "Species\tComplete_percent\tSingle_copy_percent\tDuplicated_percent\t"
        "Fragmented_percent\tMissing_percent\tBUSCO_n\tInterpretation\n"
    )
    with args.output.open("w") as out:
        out.write(header)
        for species in args.species:
            summaries = sorted((args.busco_root / species).glob("short_summary*.txt"))
            if not summaries:
                out.write(f"{species}\tNA\tNA\tNA\tNA\tNA\tNA\tMissing\n")
                continue
            match = PATTERN.search(summaries[0].read_text())
            if not match:
                out.write(f"{species}\tNA\tNA\tNA\tNA\tNA\tNA\tParse_failed\n")
                continue
            values = match.groups()
            complete = float(values[0])
            label = "Acceptable" if complete >= 65 else "Use_with_caution" if complete >= 50 else "Low_completeness"
            out.write("\t".join((species, *values, label)) + "\n")
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
