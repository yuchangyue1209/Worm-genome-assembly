#!/usr/bin/env bash
# Environment: genome_assembly. Successful route; do not rerun full-depth failure by default.
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
need "$EXTRACT/mt_candidate_reads.fastq.gz"
fresh "$FLYE"
python "$HELPER" sample "$EXTRACT/mt_candidate_reads.fastq.gz" "$FLYE/reads.random700.fastq" --number 700 --seed 42 | tee "$FLYE/sample_stats.txt"
flye --version > "$FLYE/flye_version.txt"
flye --pacbio-hifi "$FLYE/reads.random700.fastq" --out-dir "$FLYE/assembly" --threads 16
cat "$FLYE/assembly/assembly_info.txt"
