#!/usr/bin/env bash
# Regenerates a NEW checksum manifest, without changing any biological outputs.
source "$(dirname "$0")/00_config.sh"
DIR="$ANNO/results/filtered_candidates/Scotti_filtered7_v1"
cd "$DIR"
MANIFEST="SHA256SUMS.reviewed.$(date +%Y%m%d_%H%M%S)"
[[ ! -e "$MANIFEST" ]]
find . -maxdepth 1 -type f ! -name 'SHA256SUMS*' -print0 \
 | sort -z | xargs -0 sha256sum > "$MANIFEST"
sha256sum -c "$MANIFEST"
# Functional evidence has its own SHA256SUMS in functional_evidence_v1.
