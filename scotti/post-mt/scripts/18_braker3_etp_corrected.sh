#!/usr/bin/env bash
# Environment: braker3_annotation. Completed historical launch. Do not rerun an existing/completed training job.
# This script does not install packages or remove historical failed runs.
source "$(dirname "$0")/00_config.sh"
RUN="$ANNO/results/braker3_etp_v1_retry2"
export AUGUSTUS_CONFIG_PATH="$ANNO/config/augustus_scotti_v1"
export GENEMARK_PATH=/work/cyu/annotation/software/GeneMark-ETP/bin
export PROTHINT_PATH="$GENEMARK_PATH/gmes/ProtHint/bin"
# Keep BRAKER's Python first. Prepending stringtie_env can select Python without Bio.
export PATH="/home/cyu/.conda/envs/braker3_annotation/bin:/home/cyu/.conda/envs/stringtie_env/bin:$GENEMARK_PATH:$PROTHINT_PATH:/work/cyu/annotation/software/GeneMark-ETP/tools:$PATH"
hash -r
command -v python3
python3 -c 'from Bio.Seq import Seq; import Bio; print("Biopython PASS:", Bio.__version__)'
for tool in braker.pl augustus bam2hints samtools stringtie gmetp.pl prothint.py; do command -v "$tool"; done
GENOME="$ROOT/repeats/02_repeatmasker/Scotti_nuclear.mt_filtered.ragtag.provisional.fa.masked"
BAM="$ANNO/results/star/02_rna_alignment/Scotti_nuclear.rnaseq.bam"
PROTEINS="$ANNO/results/braker3_etp_v1/input.proteins.fa"
for file in "$GENOME" "$BAM" "$PROTEINS"; do [[ -s "$file" ]]; done
[[ -w "$AUGUSTUS_CONFIG_PATH/species" ]]
[[ ! -e "$AUGUSTUS_CONFIG_PATH/species/scotti_cyu_etp_v1" ]] || {
 echo 'Training species exists; inspect before restarting'; exit 1;
}
samtools quickcheck -v "$BAM"
mkdir "$RUN"
sha256sum "$GENOME" "$BAM" "$PROTEINS" > "$RUN/inputs.sha256"
braker.pl --genome="$GENOME" --bam="$BAM" --prot_seq="$PROTEINS" \
 --softmasking --species=scotti_cyu_etp_v1 --workingdir="$RUN" --threads=32 \
 --GENEMARK_PATH="$GENEMARK_PATH" --PROTHINT_PATH="$PROTHINT_PATH" --gff3 \
 2>&1 | tee "$RUN/console.log"
