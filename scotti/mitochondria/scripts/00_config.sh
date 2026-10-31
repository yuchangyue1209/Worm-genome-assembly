#!/usr/bin/env bash
# Source from stage scripts. Run on stickleback, not the local Mac.
set -euo pipefail
ROOT=${ROOT:-/work/cyu/assembly/worm/scotti}
MT="$ROOT/mitochondria"
RAW_ASM="$ROOT/raw/scotti_2_ragtag_500min_direct.fasta.scaffolds.simpleheaders.fa"
BAM="$ROOT/reads/hifi/scotti_2_merged.bam"
READS="$ROOT/reads/hifi/scotti_2.from_bam.fastq.gz"
CAND="$ROOT/excluded/scotti_mt_candidates.fa"
EXTRACT="$MT/01_read_extraction_from_bam"
FLYE="$MT/03_flye_hifi_random700"
ASM="$FLYE/assembly/assembly.fasta"
VAL="$MT/04_validation"
MITOS="$MT/05_mitos2"
SOL=/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa
SOL_GFF=/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.gff3
SOL_REF=/work/cyu/assembly/worm/mitochondria/Ssolidus_mt_ref.fa
DBBASE=/work/cyu/assembly/worm/mitochondria/mitos_db
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
HELPER="$SCRIPT_DIR/mt_helpers.py"
need() { for input in "$@"; do [[ -s "$input" ]] || { echo "Missing/empty: $input" >&2; exit 1; }; done; }
# A fresh directory prevents accidental reruns over prior results.
fresh() { mkdir -p "$(dirname "$1")"; mkdir "$1"; }
