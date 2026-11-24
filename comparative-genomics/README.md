# *Schistocephalus* Comparative Genomics

This directory contains the downstream multi-species comparative-genomics
workflow. It is intentionally separate from `../solidus/workflow-scripts/`, which
documents assembly and annotation of the *S. solidus* reference genome.

## Current species set

- `S_solidus`
- `L_intestinalis`
- `Spirometra_Aus1`
- `D_latus`
- `S_proliferum`
- `H_microstoma`
- `E_multilocularis`
- `E_granulosus`
- `T_solium`
- `T_multiceps`

These ten species remain the default background array for existing jobs.
The *S. cotti* scientific release is now published at
`/work/shared/Scotti_release_v1.0_preNCBI_20261118`, with 9,234 representative
proteins selected as the baseline. Eleven-species OrthoFinder and independent
phylogenetic analyses have already been completed with the filtered7 baseline;
see [current status and exact counts](CURRENT_STATUS.md). New script-03 runs
can explicitly include the shared release with `INCLUDE_SCOTTI_RELEASE=1`.
Existing frozen runs and standardized links remain unchanged. The fragmented
*S. erinaceieuropaei* PRJEB1202 annotation is deliberately excluded from the
primary analysis; *Spirometra* sp. Aus1 is used as the representative of that
lineage.

## Layout

```text
comparative-genomics/
  config.example.sh
  scripts/
    00_standardize_new_cyclophyllidea.sh
    00b_standardize_s_proliferum.sh
    00c_audit_s_proliferum.py
    01_input_and_protein_qc.sh
    02_busco_genomes.sh
    02b_busco_proteins.sh
    02c_busco_summary.py
    03_orthofinder.sh
    04_interproscan.sh
    05_eggnog_mapper.sh
    06_functional_coverage.sh
    07_extreme_family_qc.py
    08_deeptmhmm2.sh
    09_integrate_secretome.py
    09_integrate_secretome.sh
    10_parse_merops.py
    10_merops.sh
    11_plot_panel_a_exploratory.py
    12_plot_iqtree_concordance.py
    13_prepare_singlecopy_phylogeny.py
    14_run_iqtree_concordance.sh
    15_summarize_deep_topology_support.py
    16_run_scotti_without32_sensitivity_orthofinder.sh
```

Copy the configuration and edit paths for the server:

```bash
cp comparative-genomics/config.example.sh \
   comparative-genomics/config.sh
```

Scripts use the standardized project layout:

```text
02_standardized_data/
  genomes/<species>.fasta
  gff/<species>.gff3
  cds/<species>.fna
  proteins_primary/<species>.faa
```

## Recommended order

1. Standardize newly downloaded taxa and validate all inputs. For
   *S. proliferum*, run `00b_standardize_s_proliferum.sh` followed by
   `00c_audit_s_proliferum.py`.
2. Calculate assembly/protein QC and run genome- and protein-mode BUSCO.
3. Run OrthoFinder and infer the species tree.
4. Run uniform InterProScan and eggNOG-mapper annotation.
5. Audit extreme gene families for TE and fragmented-gene artefacts.
6. Run DeepTMHMM2 and SignalP 6, then integrate their predictions with the
   reproducible classical-secretome rule below.
7. Audit TE-associated/extreme families and document any exclusion policy.
   Recover the published *S. proliferum* TE IDs if using that comparison;
   whole-proteome TE removal is not a universal OrthoFinder requirement.
8. Use the released 9,234-record *S. cotti* baseline for new input freezes.
   Compare checksums with completed inputs before rerunning identical analyses.
   Validate the chosen family/locus inputs before CAFE5 or molecular evolution.

CAFE5 should not be run directly on the unfiltered OrthoFinder count matrix.
TE-derived families, fragmented gene models, extreme families, and poorly
annotated taxa must first be assessed.

## Environments

- OrthoFinder, IQ-TREE, BUSCO and sequence utilities: project environment
- eggNOG-mapper 2.1.13: Python-compatible eggNOG environment
- DeepTMHMM2: Python 3.11 environment (`dtm2` command)
- SignalP 6: separate Python 3.10 environment and licensed model package

Expensive searches are skipped when a validated final output already exists.

After SignalP 6 and DeepTMHMM2 are complete for every configured species, run:

```bash
bash comparative-genomics/scripts/09_integrate_secretome.sh
```

The integration uses tab-delimited SignalP output and normalizes FASTA headers
to their first token. A classical secreted protein must be SignalP-positive,
must not be classified as a beta-barrel protein, and must have no
DeepTMHMM2-predicted transmembrane helix beginning after residue 40. The cutoff
can be changed explicitly with `TM_START_CUTOFF`, but 40 is the documented
project default. Final tables, ID lists, and FASTA files are written to
`06_secretome/integrated_final/`.

Run the uniform MEROPS search and classification of the final classical
secretomes with:

```bash
bash comparative-genomics/scripts/10_merops.sh
```

The filter retains hits with E-value at most `1e-5` and subject coverage of at
least 50%, then selects the highest-bitscore retained hit per query. No 30%
identity cutoff is imposed because the validated legacy results contain
lower-identity, coverage-supported matches. Existing non-empty raw BLAST
results are reused only within the dedicated `raw_classical_secretome/`
directory. The primary-proteome search is not used for the secreted
host-interaction comparison. Combined outputs use a
`secretome_*_all_species.tsv` suffix so legacy tables and any exploratory
whole-proteome searches are not overwritten.

Combine the configured protein-mode BUSCO summaries with:

```bash
source comparative-genomics/config.sh
python comparative-genomics/scripts/02c_busco_summary.py \
  --busco-root "${QC}/busco_protein" \
  --species "${SPECIES[@]}" \
  --output "${QC}/all_species_protein_busco.tsv"
```

## Current completion boundary

The cotti v1.0 baseline is released: **9,234 genes, 11,464 protein isoforms,
9,234 longest representatives**. It retains all 32 strong TE-associated
candidate genes and the 47 unresolved cases. The 9,202-record without32 set
is sensitivity provenance only. Its 672 BUSCO group statuses match the
baseline exactly. [Release files and checks](../scotti/release_v1.0.md).

The completed Cyclophyllidea-rooted OrthoFinder run has **13,937 final
orthogroups and 138,538 assigned proteins** in mutually consistent final
membership/count tables. Its cached 374-group single-copy list is inconsistent
with those final group IDs; 373 strict one-copy-per-species groups were
reselected from final membership and matched to their alignments.

The independent 373-locus supermatrix has 190,437 amino-acid sites. Deep
gene-tree support remains conflicted: 144/373 Sparganum/Spirometra basal,
68/373 Schistocephalus basal, 129/373 D. latus/Ligula basal and 32/373 other
discordances. The cotti/solidus relationship is stable (gCF 86.06%, sCF 92.20%).
See CURRENT_STATUS.md for run-specific summary/membership discrepancies.

Release publication did not replace comparative inputs or recompute analyses.
Use a checksum comparison before rerunning completed 9,234-protein baselines.
The optional without32 analysis remains a distinct comparison in script 16;
it is not a required next step or the selected release proteome. If actual
protein inputs change, generate new OrthoFinder-dependent outputs in a new
run. Existing functional searches can be reused through protein-ID mappings.

To explicitly freeze the released baseline in a new OrthoFinder run, after
checking whether a new run is needed:

```bash
INCLUDE_SCOTTI_RELEASE=1 INCLUDE_SCOTTI_CANDIDATE=0 \
ORTHO_RUN_ID=release11_scotti_v1_new_run \
bash comparative-genomics/scripts/03_orthofinder.sh
```

This is an opt-in command, not a job submitted by this update. The default
background array for other functional/secretome scripts remains unchanged.

The first server sensitivity attempt named
`exploratory11_scotti_without32_v1` is invalid and must not be interpreted:
the staging function expanded its target before assigning the local species
name, so the cotti copy overwrote `T_multiceps.faa`; OrthoFinder consequently
received ten files. Script 16 is corrected, verifies all staged counts and
SHA256 values after copying, and defaults to the new `_v2` run directory.

Do not transfer internal-branch duplication counts from the OrthoFinder tree
onto a different IQ-TREE topology. Terminal duplication counts may be shown
separately with their exploratory status stated. Final CAFE5, focal-species
gain/loss, dN/dS, branch-site, TE–gene association, and microsynteny analyses
still require validated analysis-specific inputs; a rerun is needed when effective inputs or methods change.

The public *S. proliferum* SprJ_4.0 annotation contains 29,231 gene/mRNA
models and 28,993 translated proteins. Each gene has at most one mRNA, so no
longest-isoform selection is required. The 238 untranslated models are retained
in an audit table but are absent from protein-based analyses. The paper reports
6,316 TE-associated loci; these must be recovered from the official
supplementary data rather than inferred by keyword matching the GFF3.

Low-completeness annotations, especially *L. intestinalis*, *D. latus*, and
*T. multiceps*, must not be used as sole evidence for gene loss or family
contraction.
