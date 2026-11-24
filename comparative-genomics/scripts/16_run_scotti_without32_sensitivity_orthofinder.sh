#!/usr/bin/env bash
# OPTIONAL paired OrthoFinder sensitivity run; the selected release retains 9,234 proteins.
# This is not a prerequisite for release and is not launched by the documentation update.
# This does not overwrite the 9,234-protein exploratory run or promote the
# sensitivity FASTA to a final annotation release.
set -euo pipefail

ROOT=${ROOT:-/work/cyu/Schistocephalus_comparative_genomics}
RUN=${RUN:-"${ROOT}/04_orthofinder/runs/exploratory11_scotti_without32_v2"}
SCOTTI_RELEASE=${SCOTTI_RELEASE:-/work/shared/Scotti_release_v1.0_preNCBI_20261118}
COTTI=${COTTI:-"${SCOTTI_RELEASE}/working_snapshot/TE_sensitivity/Scotti_primary.without32_strong_TE_candidates.sensitivity.faa"}
THREADS=${THREADS:-32}
ANALYSIS_THREADS=${ANALYSIS_THREADS:-16}

command -v orthofinder >/dev/null
[[ -s "$COTTI" ]] || { echo "ERROR: missing cotti sensitivity FASTA: $COTTI" >&2; exit 1; }
[[ ! -e "$RUN" ]] || { echo "ERROR: run directory already exists: $RUN" >&2; exit 1; }

species=(
  S_solidus L_intestinalis Spirometra_Aus1 D_latus S_proliferum
  E_multilocularis H_microstoma E_granulosus T_solium T_multiceps
)

mkdir -p "$RUN/input"
printf 'species\tproteins\tsource\tsha256\n' > "$RUN/input_manifest.tsv"

stage() {
  local sp source target count digest duplicates
  sp=$1
  source=$2
  target="$RUN/input/${sp}.faa"
  [[ -s "$source" ]] || { echo "ERROR: missing input: $source" >&2; exit 1; }
  cp -p "$source" "$target"
  count=$(awk '/^>/{n++} END{print n+0}' "$target")
  duplicates=$(awk '/^>/{id=$1;sub(/^>/,"",id);seen[id]++} END{for(id in seen)if(seen[id]>1)n++;print n+0}' "$target")
  [[ "$count" -gt 0 && "$duplicates" -eq 0 ]] || { echo "ERROR: invalid FASTA for $sp" >&2; exit 1; }
  digest=$(sha256sum "$target" | awk '{print $1}')
  printf '%s\t%s\t%s\t%s\n' "$sp" "$count" "$source" "$digest" >> "$RUN/input_manifest.tsv"
}

for sp in "${species[@]}"; do
  stage "$sp" "$ROOT/02_standardized_data/proteins_primary/${sp}.faa"
done
stage S_cotti "$COTTI"

cotti_n=$(awk -F '\t' '$1=="S_cotti"{print $2}' "$RUN/input_manifest.tsv")
[[ "$cotti_n" -eq 9202 ]] || { echo "ERROR: expected 9202 cotti proteins; found $cotti_n" >&2; exit 1; }

input_n=$(find "$RUN/input" -maxdepth 1 -type f -name '*.faa' | wc -l)
[[ "$input_n" -eq 11 ]] || { echo "ERROR: expected 11 staged FASTA files; found $input_n" >&2; exit 1; }

while IFS=$'\t' read -r sp expected source expected_sha; do
  [[ "$sp" != species ]] || continue
  target="$RUN/input/${sp}.faa"
  [[ -s "$target" ]] || { echo "ERROR: staged file missing: $target" >&2; exit 1; }
  observed=$(awk '/^>/{n++} END{print n+0}' "$target")
  observed_sha=$(sha256sum "$target" | awk '{print $1}')
  [[ "$observed" -eq "$expected" ]] || { echo "ERROR: staged count mismatch for $sp" >&2; exit 1; }
  [[ "$observed_sha" == "$expected_sha" ]] || { echo "ERROR: staged hash mismatch for $sp" >&2; exit 1; }
done < "$RUN/input_manifest.tsv"

cat > "$RUN/STATUS.txt" <<'EOF'
Exploratory eleven-species TE-exclusion sensitivity analysis.
S. cotti input omits 32 strong TE-associated primary-protein candidates.
The original 9,234-protein run remains the baseline and is not overwritten.
This is not a final non-TE annotation or frozen comparative input.
EOF

orthofinder -f "$RUN/input" -t "$THREADS" -a "$ANALYSIS_THREADS" \
  2>&1 | tee "$RUN/orthofinder.log"

printf 'PASS\t%s\n' "$(date -Is)" > "$RUN/ORTHOFINDER_COMPLETE.tsv"
