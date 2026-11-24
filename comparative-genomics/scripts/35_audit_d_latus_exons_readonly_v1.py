#!/usr/bin/env python3
"""Read-only D. latus exon audit; reuses the selected Table 1 collector."""
import csv
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

sys.dont_write_bytecode = True
root = Path(sys.argv[1] if len(sys.argv) > 1 else
            "/work/cyu/Schistocephalus_comparative_genomics")
script = Path(sys.argv[2]) if len(sys.argv) > 2 else root / "comparative-genomics/scripts/30_collect_release_table1.py"
genome = root / "02_standardized_data/genomes/D_latus.fasta"
gff = root / "02_standardized_data/gff/D_latus.gff3"
protein = root / "02_standardized_data/proteins_primary/D_latus.faa"
table = root / "03_qc/release_table1_11species_v4/Table1_11species.tsv"
for path in (script, genome, gff, protein, table):
    if not path.is_file():
        raise SystemExit(f"Missing required file: {path}")

spec = importlib.util.spec_from_file_location("collector", script)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
print("=== Inputs ===", flush=True)
for path in (genome, gff, protein, table):
    print(path, "->", path.resolve(), flush=True)
print("=== Recompute statistics and compare with v4; scanning genome FASTA ===", flush=True)
assembly = m.fasta_stats(genome, nucleotide=True)
primary = m.fasta_stats(protein)
stats, problems = m.annotation_stats(gff, assembly["lengths"], set(primary["lengths"]))
with table.open() as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)
    index = header.index("D_latus")
    saved = {row[0]: row[index] for row in reader if len(row) == len(header)}
print("Metric\tRecomputed\tTable1_v4\tMatches_display_precision")
legacy_labels = {
    "exons": "Exon feature records (all annotated transcripts)",
    "coding_exons": "Exon feature records (protein-coding transcripts)",
    "mean_exon_bp": "Average protein-coding exon feature length (bp)",
}
for key, label, decimals in m.METRICS:
    if key not in stats:
        continue
    value = stats[key]
    shown = f"{value:.{decimals}f}" if value is not None else "n.a."
    old = saved.get(label, saved.get(legacy_labels.get(key), "NOT_FOUND"))
    print(label, shown, old, shown == old, sep="\t")
print("Representative_proteins", primary["records"], sep="\t")
print("Representative_mapping_problems", len(problems), sep="\t")
for item in problems[:10]:
    print("MAPPING_PROBLEM", item, sep="\t")

features = Counter()
exact_exons = Counter()
parent_exons = Counter()
raw_by_tx = Counter()
unique_by_tx = defaultdict(set)
gene_to_tx = defaultdict(set)
transcripts, cds_parents = set(), set()
with m.opengz(gff) as handle:
    for line in handle:
        if line.startswith("##FASTA"):
            break
        if not line.strip() or line.startswith("#"):
            continue
        c = line.rstrip("\r\n").split("\t")
        if len(c) != 9:
            raise SystemExit("Malformed GFF row")
        features[c[2]] += 1
        a = m.attributes(c[8])
        parents = tuple(sorted(unquote(x) for x in a.get("Parent", "").split(",") if x))
        fid = unquote(a.get("ID", ""))
        if c[2] in ("mRNA", "transcript"):
            transcripts.add(fid)
            for parent in parents:
                gene_to_tx[parent].add(fid)
        if c[2] == "CDS":
            cds_parents.update(parents)
        if c[2] == "exon":
            coord = (c[0], int(c[3]), int(c[4]), c[6])
            exact_exons[tuple(c)] += 1
            parent_exons[(coord, parents)] += 1
            for parent in parents:
                raw_by_tx[parent] += 1
                unique_by_tx[parent].add(coord)
print("=== Feature records and repeated exon records ===")
for kind, count in sorted(features.items()):
    print(kind, count, sep="\t")
print("Unique_mRNA_transcript_IDs", len(transcripts), sep="\t")
print("CDS_linked_mRNA_transcript_IDs", len(transcripts & cds_parents), sep="\t")
print("Gene_parent_IDs_with_multiple_transcripts", sum(len(v) > 1 for v in gene_to_tx.values()), sep="\t")
print("Extra_identical_exon_rows", sum(n - 1 for n in exact_exons.values()), sep="\t")
print("Extra_exon_rows_same_coordinates_and_Parent", sum(n - 1 for n in parent_exons.values()), sep="\t")
print("Transcript_ID\tRaw_exon_rows\tDistinct_exon_coordinates\tExtra_rows")
different = sorted((n - len(unique_by_tx[t]), t) for t, n in raw_by_tx.items() if n != len(unique_by_tx[t]))
print("Transcripts_with_repeated_exon_coordinates", len(different), sep="\t")
for extra, tid in sorted(different, key=lambda x: (-x[0], x[1]))[:10]:
    print(tid, raw_by_tx[tid], len(unique_by_tx[tid]), extra, sep="\t")
print("NOTE: repeated records are reported, not deleted; no source files were modified.")
