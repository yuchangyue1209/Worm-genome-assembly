#!/usr/bin/env bash
# Environment: repeat_annotation. Completed historically; do not rerun by default.
source "$(dirname "$0")/00_config.sh"
LIB="$ROOT/repeats/01_repeatmodeler_fixed/Scotti_ragtag_nuclear-families.fa"
OUT="${OUT:-$ROOT/repeats/02_repeatmasker}"
[[ -s "$NUCLEAR" && -s "$LIB" ]]
mkdir "$OUT"
sha256sum "$NUCLEAR" "$LIB" > "$OUT/inputs.sha256"
RepeatMasker -pa 8 -lib "$LIB" -xsmall -gff -dir "$OUT" "$NUCLEAR" 2>&1 | tee "$OUT/run.log"
