#!/usr/bin/env bash
# Verify that the selected S. cotti release, frozen snapshot and existing
# OrthoFinder baseline use the same 9,234 representative proteins.
set -euo pipefail

EXPECTED_PROTEINS=${EXPECTED_PROTEINS:-9234}
EXPECTED_SHA256=${EXPECTED_SHA256:-07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371}

RELEASE_FASTA=${RELEASE_FASTA:-/work/shared/Scotti_release_v1.0_preNCBI_20261118/proteins/Scotti_proteins.longest_isoform.v1.0.fa}
SNAPSHOT_FASTA=${SNAPSHOT_FASTA:-/work/shared/Scotti_release_v1.0_preNCBI_20261118/working_snapshot/protein_annotation/Scotti_nuclear.filtered7.longest.proteins.fa}
ORTHOFINDER_FASTA=${ORTHOFINDER_FASTA:-/work/cyu/Schistocephalus_comparative_genomics/04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/S_cotti.faa}

printf 'Artifact\tProteins\tSHA256\tPath\n'
for entry in \
    "Scientific_release|$RELEASE_FASTA" \
    "Frozen_working_snapshot|$SNAPSHOT_FASTA" \
    "OrthoFinder_baseline|$ORTHOFINDER_FASTA"; do
    label=${entry%%|*}
    fasta=${entry#*|}
    [[ -s "$fasta" ]] || { echo "ERROR: missing or empty $fasta" >&2; exit 1; }

    proteins=$(awk '/^>/{n++} END{print n+0}' "$fasta")
    digest=$(sha256sum "$fasta" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\n' "$label" "$proteins" "$digest" "$fasta"

    [[ "$proteins" -eq "$EXPECTED_PROTEINS" ]] || {
        echo "ERROR: $label has $proteins proteins; expected $EXPECTED_PROTEINS" >&2
        exit 1
    }
    [[ "$digest" == "$EXPECTED_SHA256" ]] || {
        echo "ERROR: $label SHA256 is $digest; expected $EXPECTED_SHA256" >&2
        exit 1
    }
done

echo "PASS: release, snapshot and OrthoFinder baseline are byte-identical."
