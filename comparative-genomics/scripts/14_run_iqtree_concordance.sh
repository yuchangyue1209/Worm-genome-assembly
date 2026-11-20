#!/usr/bin/env bash
# Run partitioned concatenation, site concordance, and gene concordance.
set -euo pipefail

: "${PHYLO_DIR:?Set PHYLO_DIR to the output of 13_prepare_singlecopy_phylogeny.py}"
THREADS=${THREADS:-20}
SEED=${SEED:-20261109}
IQTREE=${IQTREE:-iqtree3}

command -v "$IQTREE" >/dev/null
cd "$PHYLO_DIR"
[[ -s supermatrix.fa && -s partitions.nex && -d trimmed ]]

"$IQTREE" -s supermatrix.fa -p partitions.nex -st AA -m MFP \
  -B 1000 -alrt 1000 -bnni -T "$THREADS" -seed "$SEED" \
  --prefix concatenated

mkdir -p site_concordance_v1
"$IQTREE" -te concatenated.treefile -s supermatrix.fa \
  -p concatenated.best_scheme.nex --scfl 100 -T "$THREADS" -seed "$SEED" \
  --prefix site_concordance_v1/concordance

mkdir -p gene_concordance_v1
"$IQTREE" -S trimmed -st AA -m MFP -T "$THREADS" -seed "$SEED" \
  --prefix gene_concordance_v1/loci
"$IQTREE" -t concatenated.treefile \
  --gcf gene_concordance_v1/loci.treefile \
  --prefix gene_concordance_v1/concordance

printf 'PASS\t%s\n' "$(date -Is)" > IQTREE_CONCORDANCE_COMPLETE.tsv
