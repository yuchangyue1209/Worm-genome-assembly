#!/usr/bin/env bash
# Historical setup already completed. Do NOT rerun over current configuration.
# Run only to reconstruct the initial setup in a clean project.
source "$(dirname "$0")/00_config.sh"
SOURCE_CONFIG=/home/cyu/.conda/envs/braker3_annotation/config
CONFIG="$ANNO/config/augustus_scotti_v1"
RUN="$ANNO/results/braker3_etp_v1"
SOURCE_PROT=/work/cyu/annotation/ssol-annotation/results/evidence/braker_proteins.filtered.fa
[[ -d "$SOURCE_CONFIG/species" && -s "$SOURCE_PROT" ]]
[[ ! -e "$CONFIG" && ! -e "$RUN" ]]
mkdir -p "$ANNO/config" "$ANNO/results"
mkdir "$CONFIG"; cp -a "$SOURCE_CONFIG"/. "$CONFIG"/
mkdir "$RUN"; cp "$SOURCE_PROT" "$RUN/input.proteins.fa"
# Independent config and species preserve the solidus training results.
