#!/usr/bin/env bash
# Screen raw focal species-specific proteins against the sister proteome and,
# for proteins without >=50% query coverage, against the sister genome.
set -euo pipefail

ROOT=${ROOT:-/work/cyu/Schistocephalus_comparative_genomics}
RUN=${RUN:-"${ROOT}/04_orthofinder/runs/exploratory11_scotti_filtered7_v1"}
CANDIDATES=${CANDIDATES:-"${ROOT}/07_gene_family/species_specific_baseline9234"}
OUT=${OUT:-"${CANDIDATES}/reciprocal_search_v1"}
COTTI_GENOME=${COTTI_GENOME:-/work/cyu/assembly/worm/scotti/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa}
SOLIDUS_GENOME=${SOLIDUS_GENOME:-"${ROOT}/02_standardized_data/genomes/S_solidus.fasta"}
THREADS=${THREADS:-16}
EVALUE=${EVALUE:-1e-5}
QCOV=${QCOV:-50}
FORCE=${FORCE:-0}

for command in makeblastdb blastp tblastn python3; do
    command -v "$command" >/dev/null 2>&1 || {
        echo "ERROR: command not found: $command" >&2
        exit 1
    }
done

COTTI_PROTEOME="${RUN}/input/S_cotti.faa"
SOLIDUS_PROTEOME="${RUN}/input/S_solidus.faa"
COTTI_QUERY="${CANDIDATES}/S_cotti.species_specific.faa"
SOLIDUS_QUERY="${CANDIDATES}/S_solidus.species_specific.faa"

for file in \
    "$COTTI_PROTEOME" "$SOLIDUS_PROTEOME" \
    "$COTTI_QUERY" "$SOLIDUS_QUERY" \
    "$COTTI_GENOME" "$SOLIDUS_GENOME"; do
    [[ -s "$file" ]] || { echo "ERROR: missing or empty $file" >&2; exit 1; }
done

if [[ -e "$OUT" && "$FORCE" != 1 ]]; then
    echo "ERROR: output exists: $OUT (use a new OUT or FORCE=1)" >&2
    exit 1
fi
if [[ -e "$OUT" && "$FORCE" == 1 ]]; then
    stamp=$(date +%Y%m%d_%H%M%S)
    mv "$OUT" "${OUT}.previous_${stamp}"
fi

mkdir -p "$OUT/db" "$OUT/blastp" "$OUT/tblastn"

makeblastdb -in "$SOLIDUS_PROTEOME" -dbtype prot -parse_seqids \
    -out "$OUT/db/S_solidus" >/dev/null
makeblastdb -in "$COTTI_PROTEOME" -dbtype prot -parse_seqids \
    -out "$OUT/db/S_cotti" >/dev/null

blastp -query "$COTTI_QUERY" -db "$OUT/db/S_solidus" \
    -evalue "$EVALUE" -seg yes -max_target_seqs 50 -num_threads "$THREADS" \
    -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' \
    -out "$OUT/blastp/S_cotti_vs_S_solidus.tsv"

blastp -query "$SOLIDUS_QUERY" -db "$OUT/db/S_cotti" \
    -evalue "$EVALUE" -seg yes -max_target_seqs 50 -num_threads "$THREADS" \
    -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' \
    -out "$OUT/blastp/S_solidus_vs_S_cotti.tsv"

python3 - "$OUT" "$COTTI_QUERY" "$SOLIDUS_QUERY" "$QCOV" <<'PY'
from pathlib import Path
import sys

out = Path(sys.argv[1])
cotti_query = Path(sys.argv[2])
solidus_query = Path(sys.argv[3])
qcov_cutoff = float(sys.argv[4])

def ids(path):
    with path.open() as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]

def best_hits(path):
    best = {}
    if path.stat().st_size:
        with path.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                key = fields[0]
                rank = (float(fields[8]), -float(fields[7]))
                if key not in best or rank > best[key][0]:
                    best[key] = (rank, line.rstrip())
    return best

def write_pair(query_species, target_species, query_fasta):
    raw = out / "blastp" / f"{query_species}_vs_{target_species}.tsv"
    best = best_hits(raw)
    best_file = out / "blastp" / f"{query_species}_vs_{target_species}.best.tsv"
    with best_file.open("w") as handle:
        for protein_id in sorted(best):
            handle.write(best[protein_id][1] + "\n")
    unresolved = [
        protein_id
        for protein_id in ids(query_fasta)
        if protein_id not in best or float(best[protein_id][1].split("\t")[6]) < qcov_cutoff
    ]
    id_file = out / "blastp" / f"{query_species}.no_qcov{int(qcov_cutoff)}.ids"
    id_file.write_text("".join(f"{protein_id}\n" for protein_id in unresolved))
    return unresolved

def subset_fasta(source, wanted, output):
    wanted = set(wanted)
    keep = False
    with source.open() as inp, output.open("w") as out_handle:
        for line in inp:
            if line.startswith(">"):
                keep = line[1:].split()[0] in wanted
            if keep:
                out_handle.write(line)

cotti = write_pair("S_cotti", "S_solidus", cotti_query)
solidus = write_pair("S_solidus", "S_cotti", solidus_query)
subset_fasta(cotti_query, cotti, out / "blastp" / "S_cotti.no_qcov50.faa")
subset_fasta(solidus_query, solidus, out / "blastp" / "S_solidus.no_qcov50.faa")
print(f"No protein qcov50: S_cotti={len(cotti)} S_solidus={len(solidus)}")
PY

makeblastdb -in "$SOLIDUS_GENOME" -dbtype nucl -parse_seqids \
    -out "$OUT/db/S_solidus_genome" >/dev/null
makeblastdb -in "$COTTI_GENOME" -dbtype nucl -parse_seqids \
    -out "$OUT/db/S_cotti_genome" >/dev/null

tblastn -query "$OUT/blastp/S_cotti.no_qcov50.faa" \
    -db "$OUT/db/S_solidus_genome" -evalue "$EVALUE" -seg yes \
    -max_target_seqs 50 -num_threads "$THREADS" \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send qcovs evalue bitscore' \
    -out "$OUT/tblastn/S_cotti_vs_S_solidus_genome.tsv"

tblastn -query "$OUT/blastp/S_solidus.no_qcov50.faa" \
    -db "$OUT/db/S_cotti_genome" -evalue "$EVALUE" -seg yes \
    -max_target_seqs 50 -num_threads "$THREADS" \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send qcovs evalue bitscore' \
    -out "$OUT/tblastn/S_solidus_vs_S_cotti_genome.tsv"

python3 - "$OUT" <<'PY'
from pathlib import Path
import csv
import json
import sys

out = Path(sys.argv[1])
summary = []
for query, target in (("S_cotti", "S_solidus"), ("S_solidus", "S_cotti")):
    raw = out / "tblastn" / f"{query}_vs_{target}_genome.tsv"
    best = {}
    if raw.stat().st_size:
        with raw.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                rank = (float(fields[11]), -float(fields[10]))
                if fields[0] not in best or rank > best[fields[0]][0]:
                    best[fields[0]] = (rank, line.rstrip())
    best_file = out / "tblastn" / f"{query}_vs_{target}_genome.best.tsv"
    with best_file.open("w") as handle:
        for protein_id in sorted(best):
            handle.write(best[protein_id][1] + "\n")
    q50 = sum(float(value[1].split("\t")[9]) >= 50 for value in best.values())
    q80 = sum(float(value[1].split("\t")[9]) >= 80 for value in best.values())
    summary.append((query, target, len(best), q50, q80))

with (out / "reciprocal_search_summary.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["Query_species", "Target_species", "Genome_hits", "Genome_qcov50", "Genome_qcov80"])
    writer.writerows(summary)

(out / "parameters.json").write_text(json.dumps({"qcov_threshold": 50}, indent=2) + "\n")
for row in summary:
    print("\t".join(map(str, row)))
PY

printf 'PASS\t%s\n' "$(date -Is)" > "$OUT/COMPLETE.tsv"
echo "PASS: $OUT"
