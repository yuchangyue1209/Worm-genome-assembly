#!/usr/bin/env bash
# Environment: repeat_annotation. Already started on server; do not launch duplicate job.
source "$(dirname "$0")/00_config.sh"
OUT="${OUT:-$ROOT/repeats/01_repeatmodeler_fixed}"
[[ -s "$NUCLEAR" ]]
mkdir -p "$(dirname "$OUT")"; mkdir "$OUT"; cd "$OUT"
HELP=$(RepeatModeler -help 2>&1 || true)
printf '%s\n' "$HELP" > repeatmodeler_help.txt
if [[ "$HELP" == *"-threads"* ]]; then CPU_ARGS=(-threads 32)
elif [[ "$HELP" == *"-pa"* ]]; then CPU_ARGS=(-pa 8)
else echo 'Unrecognized CPU option'; exit 1; fi
sha256sum "$NUCLEAR" > input.sha256
BuildDatabase -name Scotti_ragtag_nuclear "$NUCLEAR" 2>&1 | tee build_database.log
RepeatModeler -database Scotti_ragtag_nuclear "${CPU_ARGS[@]}" -LTRStruct 2>&1 | tee repeatmodeler.log
