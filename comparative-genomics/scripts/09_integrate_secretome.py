#!/usr/bin/env python3
"""Integrate SignalP 6 and DeepTMHMM2 predictions across species."""

from __future__ import annotations

import argparse
import csv
import json
from collections import OrderedDict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-dir", required=True, type=Path)
    parser.add_argument("--signalp-dir", required=True, type=Path)
    parser.add_argument("--deeptm-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tm-start-cutoff", type=int, default=40)
    parser.add_argument("--species", nargs="+", required=True)
    return parser.parse_args()


def read_fasta(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    identifier: str | None = None
    sequence: list[str] = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if identifier in records:
                        raise ValueError(f"Duplicate FASTA ID in {path}: {identifier}")
                    records[identifier] = "".join(sequence)
                identifier = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if identifier is not None:
        if identifier in records:
            raise ValueError(f"Duplicate FASTA ID in {path}: {identifier}")
        records[identifier] = "".join(sequence)
    return records


def read_signalp(path: Path) -> dict[str, str]:
    predictions: dict[str, str] = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if not fields or fields[0].startswith("#") or len(fields) < 2:
                continue
            identifier = fields[0].split()[0]
            if identifier in predictions:
                raise ValueError(f"Duplicate SignalP ID in {path}: {identifier}")
            predictions[identifier] = fields[1].strip()
    return predictions


def segment_values(segment: object) -> tuple[str, int | None]:
    if isinstance(segment, (list, tuple)) and segment:
        label = str(segment[0]).lower()
        try:
            start = int(segment[1]) if len(segment) > 1 else None
        except (TypeError, ValueError):
            start = None
        return label, start
    if isinstance(segment, dict):
        label = str(
            segment.get("type")
            or segment.get("label")
            or segment.get("state")
            or ""
        ).lower()
        value = segment.get("start")
        try:
            start = int(value) if value is not None else None
        except (TypeError, ValueError):
            start = None
        return label, start
    return "", None


def read_deeptm(path: Path, cutoff: int) -> dict[str, dict[str, object]]:
    data = json.loads(path.read_text())
    records = data.get("predictions", []) if isinstance(data, dict) else data
    predictions: dict[str, dict[str, object]] = {}

    for record in records:
        # predictions.json includes one metadata object without an ID.
        if not isinstance(record, dict) or not record.get("id"):
            continue
        identifier = str(record["id"]).split()[0]
        protein_type = str(record.get("type", "UNKNOWN"))
        labels: list[str] = []
        tm_starts: list[int] = []

        for segment in record.get("segments", []):
            label, start = segment_values(segment)
            labels.append(label)
            if (
                start is not None
                and ("tmhelix" in label or "transmembrane" in label)
            ):
                tm_starts.append(start)

        predictions[identifier] = {
            "type": protein_type,
            "has_signal": any("signal" in label for label in labels),
            "tm_starts": tm_starts,
            "has_post_signal_tm": any(start > cutoff for start in tm_starts),
            "is_beta_barrel": (
                "beta" in protein_type.lower()
                and "barrel" in protein_type.lower()
            ),
        }
    return predictions


def write_wrapped_fasta(
    path: Path, identifiers: list[str], sequences: OrderedDict[str, str]
) -> int:
    written = 0
    with path.open("w") as output:
        for identifier in identifiers:
            sequence = sequences.get(identifier)
            if sequence is None:
                continue
            output.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                output.write(sequence[start : start + 60] + "\n")
            written += 1
    return written


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[list[object]] = []

    for species in args.species:
        fasta_path = args.protein_dir / f"{species}.faa"
        signalp_path = args.signalp_dir / species / "prediction_results.txt"
        deeptm_path = args.deeptm_dir / species / "predictions.json"
        for path in (fasta_path, signalp_path, deeptm_path):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty input: {path}")

        sequences = read_fasta(fasta_path)
        signalp = read_signalp(signalp_path)
        deeptm = read_deeptm(deeptm_path, args.tm_start_cutoff)
        sp_ids = [identifier for identifier, value in signalp.items() if value.startswith("SP")]

        master_path = args.output_dir / f"{species}.secretome_master.tsv"
        classical_ids: list[str] = []
        post_signal_tm = beta_barrel = missing_deeptm = 0

        with master_path.open("w", newline="") as output:
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "Protein_ID",
                    "SignalP_prediction",
                    "DeepTMHMM2_type",
                    "DeepTMHMM2_signal",
                    "TM_starts",
                    f"TM_start_gt{args.tm_start_cutoff}",
                    "Beta_barrel",
                    "Classical_secretome",
                ]
            )
            for identifier in sp_ids:
                topology = deeptm.get(identifier)
                if topology is None:
                    missing_deeptm += 1
                    writer.writerow([identifier, signalp[identifier], "MISSING", 0, "", 0, 0, 0])
                    continue

                has_post_tm = bool(topology["has_post_signal_tm"])
                is_beta = bool(topology["is_beta_barrel"])
                post_signal_tm += int(has_post_tm)
                beta_barrel += int(is_beta)
                is_classical = not has_post_tm and not is_beta
                if is_classical:
                    classical_ids.append(identifier)

                writer.writerow(
                    [
                        identifier,
                        signalp[identifier],
                        topology["type"],
                        int(bool(topology["has_signal"])),
                        ",".join(map(str, topology["tm_starts"])),
                        int(has_post_tm),
                        int(is_beta),
                        int(is_classical),
                    ]
                )

        ids_path = args.output_dir / f"{species}.classical_secretome.ids"
        ids_path.write_text("".join(f"{identifier}\n" for identifier in classical_ids))
        fasta_out = args.output_dir / f"{species}.classical_secretome.faa"
        fasta_count = write_wrapped_fasta(fasta_out, classical_ids, sequences)
        if fasta_count != len(classical_ids):
            raise ValueError(
                f"{species}: {len(classical_ids)} IDs but only {fasta_count} FASTA records"
            )

        summary_rows.append(
            [
                species,
                len(sequences),
                len(signalp),
                len(sp_ids),
                len(classical_ids),
                post_signal_tm,
                beta_barrel,
                missing_deeptm,
            ]
        )

    summary_path = args.output_dir / "secretome_summary.tsv"
    with summary_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "Species",
                "Primary_proteins",
                "SignalP_records",
                "SignalP_SP",
                "Classical_secretome",
                f"SP_with_TM_start_gt{args.tm_start_cutoff}",
                "SP_beta_barrel",
                "Missing_DeepTMHMM2",
            ]
        )
        writer.writerows(summary_rows)
    print(f"PASS: {summary_path}")


if __name__ == "__main__":
    main()
