#!/usr/bin/env bash
# Source after activating the environment stated in each script.
set -euo pipefail
ROOT=/work/cyu/assembly/worm/scotti
ANNO=/work/cyu/annotation/scotti-annotation
PRIMARY="$ROOT/raw/scotti_2.asm.p_ctg.fa"
RAGTAG="$ROOT/raw/scotti_2_ragtag_500min_direct.fasta"
OLD_SAMBA="$ROOT/raw/scotti_2_ragtag_500min_direct.fasta.scaffolds.fa"
NEW_SAMBA="$ROOT/scaffolding/samba_hifi_asm_m5000/input.fa.scaffolds.fa"
HIFI="$ROOT/reads/hifi/scotti_2.from_bam.fastq.gz"
MT="$ROOT/mitochondria/08_working_release/Scotti_mt_working_v1/Scotti_mitochondrial.fa"
NUCLEAR="$ROOT/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa"
R1="$ANNO/data/rnaseq/raw/PCH03_802750_S1_L002_R1_001.fastq.gz"
R2="$ANNO/data/rnaseq/raw/PCH03_802750_S1_L002_R2_001.fastq.gz"
# Historical result directories are never reused automatically by reconstruction scripts.
