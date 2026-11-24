#!/usr/bin/env bash
# Estimate an annotation-error model, then run the primary CAFE5 model.
# TREE must be a rooted, binary, ultrametric tree in real time units.
set -Eeuo pipefail

ROOT=${ROOT:-/work/cyu/Schistocephalus_comparative_genomics}
FAMILIES=${FAMILIES:-}
TREE=${TREE:-}
OUT=${OUT:?Set OUT to the CAFE5 run directory}
K=${K:-3}
CHECK_ONLY=${CHECK_ONLY:-NO}
TIME_TREE_CONFIRMED=${TIME_TREE_CONFIRMED:-NO}
VALIDATOR=${VALIDATOR:-"${ROOT}/comparative-genomics/scripts/25_validate_cafe5_tree.py"}
ERROR_MODEL=${ERROR_MODEL:-}

validate_outputs() {
    local model=Gamma suffix result_file
    [[ "$K" != 1 ]] || model=Base
    for suffix in results family_results clade_results; do
        result_file="$OUT/02_gamma_k${K}/${model}_${suffix}.txt"
        [[ -s "$result_file" ]] || { echo "ERROR: missing or empty $result_file" >&2; return 1; }
        echo "FOUND $result_file"
    done
    # Current CAFE5 source writes count/change as .tab; older outputs
    # and README examples can use .txt. Accept either actual output.
    for suffix in change count; do
        result_file="$OUT/02_gamma_k${K}/${model}_${suffix}.tab"
        [[ -s "$result_file" ]] || result_file="$OUT/02_gamma_k${K}/${model}_${suffix}.txt"
        [[ -s "$result_file" ]] || { echo "ERROR: missing or empty ${model}_${suffix}.tab/.txt" >&2; return 1; }
        echo "FOUND $result_file"
    done
    python3 - "$OUT/02_gamma_k${K}/${model}_results.txt" "$model" <<'PY'
import math
from pathlib import Path
import re
import sys

text = Path(sys.argv[1]).read_text()
model = sys.argv[2]
def number(pattern, label):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Missing {label} in model summary")
    value = float(match.group(1))
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {label}: {value}")
    return value

likelihood = number(rf"^Model {model} Final Likelihood \(-lnL\):\s*(\S+)", "-lnL")
rate = number(r"^Lambda:\s*(\S+)", "Lambda")
if rate <= 0:
    raise ValueError("Lambda must be positive")
if model == "Gamma" and number(r"^Alpha:\s*(\S+)", "Alpha") <= 0:
    raise ValueError("Alpha must be positive")
if re.search(r"^Epsilon:", text, flags=re.MULTILINE):
    if not 0 <= number(r"^Epsilon:\s*(\S+)", "Epsilon") <= 1:
        raise ValueError("Epsilon must be between zero and one")
print(f"PASS: required output files exist; finite {model} fit (-lnL={likelihood}, Lambda={rate})")
print("NOTE: output validation does not establish convergence across independent optimization runs.")
PY
}

[[ "$K" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: K must be a positive integer" >&2; exit 1; }
if [[ "$CHECK_ONLY" == YES ]]; then
    [[ -d "$OUT" ]] || { echo "ERROR: missing run directory: $OUT" >&2; exit 1; }
    if validate_outputs; then
        if [[ -f "$OUT/STATUS.txt" && ! -e "$OUT/STATUS.before_output_check.txt" ]]; then
            cp -p "$OUT/STATUS.txt" "$OUT/STATUS.before_output_check.txt"
        fi
        printf 'PASS\t%s\tExisting outputs validated; no CAFE5 rerun\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" > "$OUT/STATUS.txt"
        echo "CAFE5 existing outputs validated: $OUT"
        exit 0
    fi
    exit 1
fi

[[ "$TIME_TREE_CONFIRMED" == "YES" ]] || {
    echo "ERROR: set TIME_TREE_CONFIRMED=YES only after confirming real time units." >&2
    echo "Do not use an IQ-TREE substitutions/site tree directly." >&2
    exit 1
}
[[ -s "$FAMILIES" ]] || { echo "ERROR: missing $FAMILIES" >&2; exit 1; }
[[ -s "$TREE" ]] || { echo "ERROR: missing $TREE" >&2; exit 1; }
[[ ! -e "$OUT" ]] || { echo "ERROR: output already exists: $OUT" >&2; exit 1; }
if [[ -n "$ERROR_MODEL" && ! -s "$ERROR_MODEL" ]]; then
    echo "ERROR: supplied error model is missing or empty: $ERROR_MODEL" >&2
    exit 1
fi

if command -v cafe5 >/dev/null 2>&1; then
    CAFE_BIN=${CAFE_BIN:-cafe5}
elif command -v cafexp >/dev/null 2>&1; then
    CAFE_BIN=${CAFE_BIN:-cafexp}
else
    echo "ERROR: neither cafe5 nor cafexp is available" >&2
    exit 1
fi

python3 "$VALIDATOR" --tree "$TREE" --families "$FAMILIES"
mkdir -p "$OUT"
trap 'rc=$?; printf "FAILED\t%s\texit=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$rc" > "$OUT/STATUS.txt"; exit "$rc"' ERR
printf 'RUNNING\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" > "$OUT/STATUS.txt"
fail_run() {
    printf 'FAILED\t%s\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" > "$OUT/STATUS.txt"
    echo "ERROR: $*" >&2
    exit 1
}
run_cafe() {
    local log_file=$1
    shift
    if "$CAFE_BIN" "$@" 2>&1 | tee "$log_file"; then
        :
    else
        fail_run "CAFE5 returned an error; inspect $log_file"
    fi
    # Some CAFE5 builds return zero even when argument parsing fails.
    if grep -Eiq 'Unrecognized parameter|Unknown option|Invalid option|(^|[[:space:]])Error:' "$log_file"; then
        fail_run "CAFE5 reported an error; inspect $log_file"
    fi
}
{
    date -u +'%Y-%m-%dT%H:%M:%SZ'
    printf 'CAFE_binary\t%s\n' "$CAFE_BIN"
    "$CAFE_BIN" --version 2>&1 || true
    sha256sum "$FAMILIES" "$TREE"
} > "$OUT/run_manifest.txt"

if [[ -z "$ERROR_MODEL" ]]; then
    echo "START error-model estimation: $(date -u)"
    run_cafe "$OUT/01_error_model.log" \
        -i "$FAMILIES" \
        -t "$TREE" \
        -o "$OUT/01_error_model" \
        -p \
        -e

    [[ -d "$OUT/01_error_model" ]] || fail_run "Missing error-model output directory"
    ERROR_MODEL=$(find "$OUT/01_error_model" -type f -name 'Base_error_model.txt' -print -quit)
    [[ -n "$ERROR_MODEL" && -s "$ERROR_MODEL" ]] || fail_run "CAFE5 did not produce Base_error_model.txt"
else
    echo "REUSE completed error model: $ERROR_MODEL"
fi
sha256sum "$ERROR_MODEL" >> "$OUT/run_manifest.txt"

echo "START gamma-k${K} analysis with estimated error model: $(date -u)"
run_cafe "$OUT/02_gamma_k${K}.log" \
    -i "$FAMILIES" \
    -t "$TREE" \
    -o "$OUT/02_gamma_k${K}" \
    -p \
    -k "$K" \
    "-e${ERROR_MODEL}"

validate_outputs || fail_run "Output validation failed; inspect the files listed above"

printf 'PASS\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" > "$OUT/STATUS.txt"
echo "CAFE5 pipeline COMPLETE: $OUT"
