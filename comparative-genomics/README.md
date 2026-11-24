# *Schistocephalus* Comparative Genomics

This directory contains the downstream multi-species comparative-genomics
workflow. It is intentionally separate from `../workflow-scripts/`, which
documents assembly and annotation of the *S. solidus* reference genome.

Start with [CURRENT_STATUS.md](CURRENT_STATUS.md) for the recorded analysis
state, selected result versions, manuscript-safe counts and unresolved checks.
This README describes how to run the workflow; having a script is not evidence
that its current eleven-species analysis has completed. Local organization
does not upload files, launch server jobs or replace frozen input data.

## Current species set

- `S_solidus`
- `S_cotti`
- `L_intestinalis`
- `Spirometra_Aus1`
- `D_latus`
- `S_proliferum`
- `H_microstoma`
- `E_multilocularis`
- `E_granulosus`
- `T_solium`
- `T_multiceps`

These eleven species form the current comparative dataset. *S. cotti* is
represented by the selected pre-NCBI scientific release
`Scotti_release_v1.0_preNCBI_20261118`, containing 9,234 representative
proteins. This protein set is frozen for the current comparative baseline;
NCBI submission and accessioning remain separate future release steps. The
fragmented
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
    00d_stage_scotti_candidate.sh
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
    17_plot_orthology_summary_clean.py
    18_plot_gene_sharing_tree_bars.py
    19_extract_focal_species_specific.py
    20_run_focal_reciprocal_searches.sh
    21_summarize_focal_species_specific.py
    22_validate_priority_candidates.py
    23_verify_scotti_release_baseline.sh
    24_prepare_cafe5_input.py
    25_validate_cafe5_tree.py
    26_run_cafe5.sh
    27_map_flagged_proteins_to_orthogroups.py
    28_make_cafe5_timetrees.R
    29_summarize_cafe5.py
    30_collect_release_table1.py
    31_review_cafe_priority_families.py
    32_plot_cafe5_tree_pies.py
    35_audit_d_latus_exons_readonly_v1.py
    36_audit_comparative_outputs_readonly.py
    audit_final_busco_inputs.py
    run_missing_genome_busco.py
    verify_missing_genome_runner.py
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

**Manuscript baseline (24 November 2026):** use the selected S. solidus 0831
and S. cotti 1118 releases, the eleven-species 9,234-cotti-protein OrthoFinder
baseline, and metazoa_odb12. Figure 1 uses 373 loci / 190,437 sites and the
16/20 focal species-specific OG counts. Figure 2 uses the existing eight-species
CAFE5 fit and all 6,628 reconstructed families; nominal family-P<0.01 screening
is a separate result, not a significance filter on Figure 2. The historical
without32 tree used for CAFE dating remains documented; do not replace it
silently with the Figure 1 tree. The steps below describe provenance and
possible future work, not an instruction to rerun completed analyses.

### Corrected entry points

- `00d_stage_scotti_candidate.sh` uses the selected shared release. Raw proteins
  are the 11,464 all-isoform set; primary proteins are the 9,234 representatives.
  It verifies representative sequences against the all-isoform set and checks
  all destinations before linking. An existing byte-identical destination is
  reused; a conflicting or broken link stops the script without overwriting it.
  Do not delete an old link just to bypass this check: audit its provenance first.
- `02e_collect_busco_summaries.py` defaults to all eleven species and odb12.
  Explicit `--lineages metazoa_odb10 metazoa_odb12` retains the historical inventory.
  It does not replace the release-aware Table 1 collector or independently
  establish that a discovered summary belongs to the selected release.
- `04_interproscan.sh` delegates to `interpro_verified_run.py`. A separate search
  FASTA removes only one terminal `*`; internal stops and malformed IDs stop the
  run. Release/OrthoFinder inputs are unchanged. Input and search hashes and
  per-protein MD5/length are recorded. New jobs require exit code zero and a
  validated TSV before a completion record is written; only hash-matching
  completed jobs are skipped. Existing historical TSVs without this record are
  preserved and cause a stop, not a rerun. The historical annotations used in the
  manuscript remain evidence; use script 31 for their focal sequence audit.
  `INTERPRO_OUT` and `INTERPRO_LOG` can select fresh directories if a new analysis
  is independently warranted. The executable hash is checked, not a complete
  fingerprint of its installed databases.

1. Standardize newly downloaded taxa and validate all inputs. For
   *S. proliferum*, run `00b_standardize_s_proliferum.sh` followed by
   `00c_audit_s_proliferum.py`.
2. Calculate assembly/protein QC and run genome- and protein-mode BUSCO.
3. Stage the selected *S. cotti* pre-NCBI release, verify its checksum against
   the frozen snapshot and OrthoFinder input, run the eleven-species
   OrthoFinder baseline, and infer the species tree.
4. Run uniform InterProScan and eggNOG-mapper annotation.
5. Audit extreme gene families for TE and fragmented-gene artefacts.
6. Run DeepTMHMM2 and SignalP 6, then integrate their predictions with the
   reproducible classical-secretome rule below.
7. Optional future sensitivity work: recover the published *S. proliferum*
   TE-associated locus list and document a non-TE protein set. This is not part
   of the existing eight-species manuscript CAFE analysis.
8. Audit focal raw species-specific orthogroups with reciprocal sister-species
   protein and genome searches, followed by RNA/repeat validation of the small
   manually prioritized set.
9. Reuse the frozen inputs and completed CAFE5 fits. Family trees, microsynteny
   and molecular-evolution analyses require separately validated inputs and
   must not be represented as completed manuscript results without evidence.

CAFE5 should not be run directly on the unfiltered OrthoFinder count matrix.
TE-derived families, fragmented gene models, extreme families, and poorly
annotated taxa must first be assessed.

## CAFE5 species sets and workflow

The recommended current screening set contains eight species:
`S_solidus`, `S_cotti`, `Spirometra_Aus1`, `E_multilocularis`,
`H_microstoma`, `E_granulosus`, `T_solium`, and `T_multiceps`.
`D_latus` and `L_intestinalis` were excluded from this screening because
their protein BUSCO odb12 completeness is only 20.4% and 28.0%, respectively;
their inclusion could confound apparent contractions. `S_proliferum` was also
omitted pending TE-family review. A nine-species run after that review and a
run excluding `T_multiceps` (57.1% protein BUSCO completeness) are possible future
sensitivity analyses, not the selected manuscript dataset or completed runs.
The complete eleven-species matrix is descriptive only.
These profiles are recorded in `data/cafe5_species_sets.tsv`.

Prepare the current eight-species screening matrix from the matching
9,234-protein OrthoFinder baseline. The following documents the completed
screen, not a synchronization step: reuse existing validated inputs instead
of rebuilding or refitting solely because code/docs were updated. The TE-OG
exclusion list must be present before preparing this matrix.

```bash
ROOT=/work/cyu/Schistocephalus_comparative_genomics
RESULT="$ROOT/04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/OrthoFinder/Results_Nov09"
PREP="$ROOT/07_gene_family/cafe5/screening8_input_v1"

python3 /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/24_prepare_cafe5_input.py \
  --gene-counts "$RESULT/Orthogroups/Orthogroups.GeneCount.tsv" \
  --out "$PREP" \
  --exclude-family-ids /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/scotti_strong32_orthogroups_v2/exclude_orthogroups.ids \
  --species \
    S_solidus S_cotti Spirometra_Aus1 \
    E_multilocularis H_microstoma E_granulosus T_solium T_multiceps
```

The script removes families found in fewer than two retained species and
separates families having at least 100 copies in any retained species. The
`families_lt100.tsv` file is used for parameter estimation. Reviewed TE or
artefact orthogroups can be removed reproducibly by adding one or more
`--exclude-family-ids path/to/IDs.txt` options. Large families are preserved
for separate manual/fixed-lambda analysis rather than silently discarded.

Keep the official 9,234-protein *S. cotti* release as the baseline, but map the
32 reviewed strong-TE protein IDs to their orthogroups and exclude/flag those
families in CAFE5 rather than altering the release FASTA:

```bash
python3 /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/27_map_flagged_proteins_to_orthogroups.py \
  --orthogroups "$RESULT/Orthogroups/Orthogroups.tsv" \
  --ids S_cotti=/work/cyu/Schistocephalus_comparative_genomics/07_gene_family/species_specific_baseline9234/S_cotti.removed32_TE.ids \
  --out "$ROOT/07_gene_family/cafe5/scotti_strong32_orthogroups_v2"
```

Then add this option to script 24:

```text
--exclude-family-ids /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/scotti_strong32_orthogroups_v2/exclude_orthogroups.ids
```

The same mapper can be used for *S. proliferum* after its official
TE-associated locus list has been converted to the exact protein IDs used by
OrthoFinder. Do not rely on loose annotation-keyword matching for that step.
The mapper automatically checks `Orthogroups_UnassignedGenes.tsv`. Flagged
proteins found there are reported separately and do not require family
exclusion because they have no row in the CAFE5 count matrix. Only IDs absent
from both the assigned and unassigned OrthoFinder tables cause a validation
error.

CAFE5 requires a rooted, binary, ultrametric **time tree** whose branch lengths
are in real time units. The IQ-TREE `concatenated.treefile` has branch lengths
in substitutions per site and must not be used directly. After producing and
documenting a calibrated time tree containing exactly the selected species,
run:

```bash
FAMILIES=/work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/screening8_input_v1/families_lt100.tsv \
TREE=/work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/time_trees/screening8_dated_16.5Ma.nwk \
OUT=/work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/runs/screening8_16.5Ma_k3_v2 \
TIME_TREE_CONFIRMED=YES \
bash /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/26_run_cafe5.sh
```

This first estimates an annotation-error model and then runs the three-category
gamma model with family p-values. Numeric ultrametricity and exact species
matching are checked before CAFE5 starts. A statistically estimated error
model is useful, but it does not rescue severely incomplete annotations and is
not a substitute for excluding the two low-completeness taxa.

The runner accepts `ERROR_MODEL=/absolute/path/Base_error_model.txt` to reuse
an estimated model with the same count matrix and time tree in an independent
fit. CAFE5 requires the filename attached to `-e`; script 26 handles this.
`CHECK_ONLY=YES OUT=/absolute/path/run K=3 bash /absolute/path/26_run_cafe5.sh`
validates existing outputs without running CAFE5. Both `.tab` and `.txt`
count/change tables are accepted. A validation pass establishes output
availability and finite parameters, not independent-run convergence.

After comparing independent fits, screen family-wide tests and describe focal
copy changes (script 29 schema v2, corrected on 2026-11-23):

```bash
python3 /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/29_summarize_cafe5.py \
  --results /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/runs/screening8_16.5Ma_k3_v2/02_gamma_k3 \
  --families /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/screening8_input_v1/families_lt100.tsv \
  --species S_cotti S_solidus \
  --family-p-threshold 0.01 \
  --fdr-threshold 0.05 \
  --pvalue-floor 0.001 \
  --out /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/summaries/screening8_16.5Ma_k3_v2_focal_nominal_p0.01_schema2_v2
```

The nominal candidate list requires a family-wide p-value below 0.01 AND a
nonzero reconstructed focal copy change. It describes a changing focal branch
within a family flagged by the tree-wide test, not a branch-specific significance
test. The threshold is an exploratory choice made after examining the first
results; do not describe it as prespecified. Re-screening at another threshold
does not require refitting the CAFE model; use a new summary directory.

Gamma branch probabilities are posterior probabilities, not p-values
([CAFE5 developer clarification](https://github.com/hahnlab/CAFE5/discussions/232)).
Schema v2 preserves these as `Branch_probability_reported`, never thresholds
them, and retires the old `Branch_pvalue_raw`, `branch_p005` and
`familyFDR_branchP` labels. Old candidate counts based on that filter must not
be cited as branch-significant results. Raw CAFE outputs and old summaries
are preserved; no OrthoFinder or CAFE fit is rerun by this summarizer.

Outputs include `focal_nominal_candidates.tsv`, `focal_family_changes.tsv`,
`focal_summary.tsv`, `family_tests.tsv`, separate `focal_BH_reported_candidates.tsv`
and `focal_BH_floor_sensitivity_candidates.tsv`, and labelled orthogroup ID
lists. Nominal p and supplementary BH q thresholds are independent options
(defaults 0.01 and 0.05). BH correction includes every finite family p-value in
the model output, not just focal or nominally selected families. Families
filtered internally by CAFE5 are reported separately and excluded from the BH
denominator. None of these family-level calculations establishes branch-level
FDR control.

CAFE5 can report zero family p-values at finite simulation/output precision.
Script 29 records their number; a reported BH q=0 is not exact-zero evidence.
The optional `--pvalue-floor 0.001` scenario applies a floor before a SEPARATE BH
calculation. Original family p-values and their BH q-values remain unchanged;
the nominal candidate list always uses the original family p-values. Without
that option the floor-sensitivity fields are NA. A floor neither recovers true
p-values nor provides a rigorous confidence bound. Keep both diagnostics;
formal FDR significance claims require valid p-values and review of model and
time-calibration assumptions.

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

Generate the clean large-font four-panel orthology figure by pairing a matching
IQ-TREE tree and OrthoFinder result directory. Do not mix the 9,234-protein
baseline tree with the 9,202-protein sensitivity statistics, or vice versa.

```bash
python3 comparative-genomics/scripts/17_plot_orthology_summary_clean.py \
  --tree /path/to/matching/concatenated.treefile \
  --results /path/to/matching/OrthoFinder/Results_* \
  --out 04_orthofinder/figures/orthology_summary_clean \
  --font-scale 1.50
```

The clean figure contains only the four panel headings, axis labels and plotted
values. It does not add an overall title, legend, explanatory subtitle,
footnote or caption file.
The revised Figure 1 uses abbreviated binomials except `Sparganum proliferum`,
which retains its full genus name at the user's request.
`Spirometra sp. Aus1` retains the genus name, with only
`Spirometra` italicized; `sp. Aus1` is roman. Default dimensions are 18.5 by
9.6 inches with larger text (`--font-scale 1.50`). Labels may use additional
blank horizontal space. Long branches over 0.20 substitutions/site are
displayed at length 0.085 and marked with double slashes; their actual and
displayed lengths are exported to `orthology_summary_clean.branch_breaks.tsv`.
The 0.05 scale applies to unbroken branch segments. Dashed tip extensions align
the names and are not evolutionary branch lengths. `--no-branch-breaks`
disables compression. The source tree and statistics are never changed.
Existing output directories are refused. A
`orthology_summary_clean.source_data.tsv` export and input hashes in
`orthology_summary_clean.audit.json` preserve numeric and typography provenance.

## Frozen *S. cotti* comparative baseline

The selected scientific release is:

```text
/work/shared/Scotti_release_v1.0_preNCBI_20261118
```

The following three 9,234-protein FASTA files were verified as byte-identical
on 2026-11-20:

```text
07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371  /work/shared/Scotti_release_v1.0_preNCBI_20261118/working_snapshot/protein_annotation/Scotti_nuclear.filtered7.longest.proteins.fa
07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371  /work/shared/Scotti_release_v1.0_preNCBI_20261118/proteins/Scotti_proteins.longest_isoform.v1.0.fa
07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371  /work/cyu/Schistocephalus_comparative_genomics/04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/S_cotti.faa
```

Consequently, publication of the shared release does not require rerunning the
baseline OrthoFinder analysis, species tree, focal species-specific analysis,
or matching descriptive figures. The historical run-directory name retains
the word `exploratory`, but its *S. cotti* protein content is the selected
scientific-release baseline. Run
`scripts/23_verify_scotti_release_baseline.sh` to repeat this identity check.
The release retains all 32 strong TE-associated candidates; their omission in
the 9,202-protein set is deliberately limited to sensitivity analysis.

## Current completion boundary

The current descriptive baseline is
`exploratory11_scotti_filtered7_v1/Results_Nov09`, containing the selected
9,234-protein *S. cotti* scientific-release baseline. It contains 16,583
orthogroups, assigns
137,611 proteins, includes 1,460 orthogroups with all eleven species, and has
373 strict one-copy-per-species orthogroups. These are reported final-membership
summaries. A saved historical Statistics_Overall.tsv instead reports 16,669
OGs / 137,591 assigned / 1,445 all-species OGs. That discrepancy has not been
reconciled; do not mix these global summaries or present them as a single
independently verified total. The clean orthology figures and
the focal species-specific analysis must use this matching baseline result and
its matching tree.

The reversible cotti TE-exclusion sensitivity input contains 9,202 proteins
after omitting 32 strong TE-associated primary candidates. The valid run is
`exploratory11_scotti_without32_v2/Results_Nov12`; it contains 16,582
orthogroups, assigns 137,750 proteins, includes 1,461 all-species orthogroups,
and has 381 strict one-copy-per-species orthogroups. Its 381-locus amino-acid
supermatrix contains 197,179 sites, and the IQ-TREE/concordance pipeline is
complete. The near-identical global summaries show that the 32-record omission
does not materially alter the broad orthology result. Never mix the baseline
tree/statistics with the sensitivity membership table, or vice versa.

The first server sensitivity attempt named
`exploratory11_scotti_without32_v1` is invalid and must not be interpreted:
the staging function expanded its target before assigning the local species
name, so the cotti copy overwrote `T_multiceps.faa`; OrthoFinder consequently
received ten files. Script 16 is corrected, verifies all staged counts and
SHA256 values after copying, and defaults to the new `_v2` run directory.

Do not transfer internal-branch duplication counts from the OrthoFinder tree
onto a different IQ-TREE topology. Terminal duplication counts may be shown
separately with their descriptive status stated. The *S. cotti* baseline is
now frozen. An eight-species CAFE5 screening and independent repeat fits are
reported complete, as recorded in CURRENT_STATUS.md; final CAFE5,
focal-species gain/loss, dN/dS, branch-site,
TE–gene association, and microsynteny analyses still require final input/QC
decisions for the remaining species and analysis-specific filtering.

The public *S. proliferum* SprJ_4.0 annotation contains 29,231 gene/mRNA
models and 28,993 translated proteins. Each gene has at most one mRNA, so no
longest-isoform selection is required. The 238 untranslated models are retained
in an audit table but are absent from protein-based analyses. The paper reports
6,316 TE-associated loci; these must be recovered from the official
supplementary data rather than inferred by keyword matching the GFF3.

Low-completeness annotations, especially *L. intestinalis*, *D. latus*, and
*T. multiceps*, must not be used as sole evidence for gene loss or family
contraction.

## Focal species-specific orthogroup QC

The raw 9,234-protein baseline contains 16 *S. cotti*-only orthogroups (74
proteins) and 20 *S. solidus*-only orthogroups (97 proteins). These are
OrthoFinder clustering summaries, not automatically novel gene families, and
the values 16 and 20 remain the correct raw values for the baseline figure.

Reciprocal sister-proteome BLASTP followed by sister-genome TBLASTN reduced the
unresolved set to seven *S. cotti* proteins and eighteen *S. solidus* proteins.
Manual model/domain review then prioritized one *S. cotti* orthogroup and four
*S. solidus* orthogroups for RNA/repeat validation. `OG0002484` is instead a
shared Fibronectin type III-containing family with a possible *S. solidus*
copy-number expansion; it is not a novel species-specific family.

Both proteins in the prioritized *S. cotti* `OG0015711` family are
transcriptionally supported (one strong and one moderate), have zero CDS
overlap with the cotti RepeatMasker annotation, and occur on the same strand
within a 94.9-kb region of scaffold 8 (83.4 kb between gene boundaries). It is
therefore a transcription-supported candidate pair, not proof of absolute
absence from every unsampled species. The scaffold is reference-guided:
whether both genes lie in one original contig has not been verified, so their
scaffold spacing is not independent tandem-duplication evidence.
In *S. solidus*, ten proteins from `OG0016342`, `OG0016337` and
`OG0009067` have strong RNA support. The three `OG0012892` models have weak or
partial support and remain provisional.

Reproduce this focal QC with versioned output directories:

```bash
BASE=/work/cyu/Schistocephalus_comparative_genomics
RESULTS="$BASE/04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/OrthoFinder/Results_Nov09"
RUN="$BASE/04_orthofinder/runs/exploratory11_scotti_filtered7_v1"
RAW="$BASE/07_gene_family/species_specific_baseline9234_v2"

python3 comparative-genomics/scripts/19_extract_focal_species_specific.py \
  --results "$RESULTS" \
  --input-dir "$RUN/input" \
  --out "$RAW" \
  --expected S_cotti=16 \
  --expected S_solidus=20

CANDIDATES="$RAW" \
OUT="$RAW/reciprocal_search_v1" \
bash comparative-genomics/scripts/20_run_focal_reciprocal_searches.sh

python3 comparative-genomics/scripts/21_summarize_focal_species_specific.py \
  --candidates "$RAW" \
  --search-dir "$RAW/reciprocal_search_v1" \
  --out "$RAW/evidence_v1" \
  --interpro S_cotti=/work/cyu/annotation/scotti-annotation/results/functional/03_interproscan/scotti_braker.tsv \
  --interpro S_solidus="$BASE/05_functional_annotation/interpro/S_solidus.tsv"

python3 comparative-genomics/scripts/22_validate_priority_candidates.py \
  --membership "$RAW/focal_species_specific_orthogroups.tsv" \
  --priority-table comparative-genomics/data/focal_species_specific_manual_review.tsv \
  --datasets comparative-genomics/data/focal_candidate_validation.server.tsv \
  --out "$RAW/priority_validation_v1"
```

An optional reviewed TE list can be added to script 21 as
`--te-ids S_cotti=/path/to/nonempty.ids`. RNA support validates transcription
but does not by itself establish lineage specificity. The evidence-screened
candidate counts belong in a supplementary validation table; they must not
silently replace the raw OrthoFinder 16/20 values in the descriptive figure.

## Eleven-species release Table 1

`scripts/30_collect_release_table1.py` collects assembly, annotation and both
genome/protein BUSCO assessments from selected inputs. It launches no new
analyses, preserves source files and refuses to overwrite an output directory.
Missing evidence stays `n.a.` with reasons in `missing_metrics.tsv`.

Exon method `transcript_exon_associations_v2` counts distinct
`(transcript ID, sequence, start, end, strand)` associations. Exon rows
attached directly to genes or pseudogenes do not count a second time when
the same interval is also attached to mRNA. Shared exons count once for each
transcript, while identical coordinates within one transcript count once.
Coding exon totals and mean lengths use explicit CDS-supported transcripts
linked to gene features and exclude pseudogene ancestry, biotypes and flags.
RNA transcript features are retained in the all-transcript total. Missing
transcript-level exons are not reconstructed from gene-level rows or CDS;
coding totals/mean lengths remain unavailable if coding transcripts lack
exon features. Per-species raw exon records, Parent levels, deduplication,
pseudogene transcript counts and missing-exon checks are recorded in
`exon_transcript_audit.tsv`; `collection.json` records the method version.

The user-supplied D. latus audit found 111,669 raw exon rows split across
gene (55,608), mRNA (55,835) and pseudogene (226) Parents. The transcript-only
audit recovered 55,835 all-transcript associations, 55,608 coding associations,
227 pseudogene-transcript associations, 2.79 exons per coding transcript and
a coding exon mean length of 227.9 bp. These are real user-supplied audit
results, not a new local/server run of all eleven genomes. Earlier v4 tables
retain their original row-count definitions and must not be overwritten.

The solidus RepeatMasker report predates removal of the 6,055-bp fish-derived
`scaffold_61`. The collector now reads that original FASTA and verifies that
the selected 0831 release differs only by the explicitly excluded sequence.
Every retained sequence must have the same ID and normalized sequence hash.
It validates all original RepeatMasker coordinates, omits `scaffold_61` hits,
and divides the retained interval union by the selected assembly length.
It does not merely copy the old rounded 68.64% from the summary report.
Exact masked/excluded bp and input fingerprints are written to
`repeatmasker_release_audit.tsv`. Failed sequence verification leaves the
repeat percentage unavailable and records the reason rather than guessing.

Run on stickleback with full paths and a new output directory:

```bash
python3 /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/30_collect_release_table1.py \
  --root /work/cyu/Schistocephalus_comparative_genomics \
  --inputs /work/cyu/Schistocephalus_comparative_genomics/03_qc/release_table1_11species_v4/selected_sources.json \
  --audit-dir /work/cyu/Schistocephalus_comparative_genomics/03_qc/final_busco_audit_20261120_194629 \
  --busco-dir /work/cyu/Schistocephalus_comparative_genomics/03_qc/final_genome_busco_odb12_v1 \
  --out /work/cyu/Schistocephalus_comparative_genomics/03_qc/release_table1_11species_v5_exonfix
```

These defaults include the solidus projection. If reusing `--inputs` from an
older collection, update its solidus `repeatmasker_out`, `repeatmasker_genome`
and `repeatmasker_exclude_ids` fields explicitly; saved configurations are
not silently overridden. TE-protein counts still require uniformly classified
representative-protein ID lists. Missing E. multilocularis exon features are
not replaced by CDS segments. Neither of these gaps triggers BUSCO, OrthoFinder,
CAFE or structural reannotation.

## Per-protein review of four CAFE candidate families

`scripts/31_review_cafe_priority_families.py` reviews `OG0000875`, `OG0000880`,
`OG0000999` and `OG0001937` from the selected 9,234-protein OrthoFinder baseline.
The reported focal membership is 26 cotti plus 11 solidus proteins. It writes
per-protein lengths, explicit release-GFF gene/transcript/CDS coordinates,
Pfam coverage and individual feature coordinates, alongside an all-eleven-species
copy-count matrix and species-prefixed per-family FASTAs for later alignment.

InterPro matches are accepted only when their sequence MD5, amino-acid length
and match coordinates agree with the frozen OrthoFinder protein input. One
terminal `*` marker may be omitted for this comparison, as verified for cotti
`g4713.t1` (81 stored characters versus 80 amino acids; the 80-aa MD5 matches
the original InterPro search input). This is recorded per protein in
`protein_sequence_match_audit.tsv`; source files, OrthoFinder inputs and
exported family FASTAs are not rewritten. Length/coverage statistics exclude
that terminal marker. Internal/multiple stops and residue substitutions are
not normalized away. A mismatched target row still stops before output creation. Missing rows
are explicitly unavailable, not evidence that a protein lacks domains.
The presence/absence of the RNR `IPR030475` active-site entry is a reported
annotation feature, not a biochemical activity prediction. Identical protein
sequences are flagged, never automatically removed. Domain architecture,
RNA support, repeat overlap and gene-tree validation remain separate checks.

```bash
python3 /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/31_review_cafe_priority_families.py \
  --root /work/cyu/Schistocephalus_comparative_genomics \
  --results /work/cyu/Schistocephalus_comparative_genomics/04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/OrthoFinder/Results_Nov09 \
  --input-dir /work/cyu/Schistocephalus_comparative_genomics/04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input \
  --out /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/candidate_review/priority4_baseline9234_v2
```

Main outputs are `family_review_summary.tsv`, `priority_protein_review.tsv`,
`feature_hits.tsv`, `family_copy_counts_allspecies.tsv`, and
`all_species_fastas/`. `family_membership.tsv`,
`priority_family_selection.tsv` and the two `*.priority_cds.bed` files prepare
later RNA/repeat work; this script does not run it or reuse results for other
candidate sets. Inputs can be overridden with `--gff SPECIES=PATH` or
`--interpro SPECIES=PATH`. Output directories are never overwritten.

This is the 37-protein CAFE priority set, not the 15-protein species-specific
candidate set validated earlier. The earlier MD5 stop was diagnosed and the
single-terminal-stop comparison fix is local and tested; successful complete
post-fix server output has not yet been supplied. A new review directory
avoids colliding with any previous partial attempt.

## CAFE time-tree figure with branch pies

`scripts/32_plot_cafe5_tree_pies.py` draws the eight-species screening result as
a time tree with large labels and red-green color-vision-friendly colors:
blue (`#0072B2`) for expanded, vermillion (`#D55E00`) for contracted, and
light gray (`#D9D9D9`) for unchanged families. The blue/vermillion pair follows
[Okabe-Ito Color Universal Design](https://jfly.uni-koeln.de/color/).
Signed `+/-` counts and named legend entries preserve meaning without color;
the two focal species use bold black labels. The palette is recorded in the
figure provenance. It writes vector PDF/SVG plus PNG,
without an overall title, long caption, invented confidence bars, WGT markers,
root copy count or pooled-change pie. All changes are shown by default, not
only nominal or FDR-screened families; the numbers are family counts, not
total gene copies gained/lost. A root pie is omitted because the root has no
incoming branch.

The script reads `Gamma_asr.tre` to recover explicit CAFE node IDs and maps
them by rooted descendant-species sets to the selected dated tree. It checks
topology, branch lengths, ultrametricity, tip labels and reconstruction
count/change/parent consistency. Available `Gamma_clade_results` counts must
agree with counts derived from positive/negative entries in `Gamma_change`.
ASR single-family copy numbers are never substituted for genome-wide totals.
No model fit or other analysis is rerun, and inputs/old figures are preserved.

```bash
/home/cyu/.conda/envs/deeptmhmm2/bin/python \
  /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/32_plot_cafe5_tree_pies.py \
  --results /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/runs/screening8_16.5Ma_k3_v2/02_gamma_k3 \
  --tree /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/time_trees/screening8_dated_16.5Ma.nwk \
  --out /work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/figures/screening8_16.5Ma_k3_v2_all_changes_cvd_v1 \
  --font-scale 1.15
```

Expected focal all-change labels are cotti `+265/-611` and solidus `+178/-230`
for the completed v2 model. These values are read from the tables, not
hard-coded. Main figure name: `cafe5_tree_pies.pdf` (also `.svg` and `.png`).
`branch_family_counts.tsv`, `node_mapping.tsv`, `figure_provenance.json` and
`interpretation_notes.txt` record definitions and exact input checksums.
Missing ASR data stops the figure rather than guessing internal node IDs;
`--asr` can select its actual path. Matplotlib is the only non-stdlib Python
dependency. Local regression/layout checks use clearly synthetic fixtures
only. The user has separately supplied the real eight-species server figure;
its raw model/figure files are not mirrored locally, and source checksums
still require retrieval from the server. Do not refit CAFE just to organize
or re-export this figure.

Optional `--family-p-threshold 0.01` produces a separately labelled nominal
family-wide subset, with pies normalized to that subset. Use a separate output
directory. This is not a branch-specific significance test or an FDR result;
branch posterior tables are never used as p-values. The supplied 16.5-Ma tree
retains its secondary-calibration and without32-tree/baseline-count provenance
limits; the axis does not establish dating accuracy.

CAFE format/interpretation references:
https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md
https://github.com/hahnlab/CAFE5/discussions/232

## Read-only synchronization and output inventory

The restored [CURRENT_STATUS.md](CURRENT_STATUS.md) records reported results,
not automatic job status. Script 36 can inspect known paths and compare small
code/document/data files against the synchronization manifest. It does not
invoke the CAFE runner's CHECK_ONLY mode (which updates STATUS.txt), source
shell configurations, modify links, or rerun searches and analyses.

```bash
python3 /work/cyu/Schistocephalus_comparative_genomics/comparative-genomics/scripts/36_audit_comparative_outputs_readonly.py \
  --root /work/cyu/Schistocephalus_comparative_genomics
```

`PRESENT_NOT_VALIDATED` means only that an expected file is nonempty. Code
hash equality is not output validation, and parsed fit parameters are not
proof of convergence. Missing paths may reflect work not done or a different
actual output location; inspect before rerunning. No live server inventory
has been performed by this local organization task.

The comparative-only archive includes source scripts, tests, configuration
examples, small input manifests and documentation. It excludes private
config.sh, Python caches, raw genomes/proteins, existing result directories
and local synthetic regression outputs. Extract to a separate directory and
compare first; do not overwrite server configuration or frozen datasets as
a consequence of a code synchronization.

The 2026-11-24 GitHub synchronization preserves lineage-separated outputs in
`scripts/02b_busco_proteins.sh` and immutable fresh-run inputs in
`scripts/03_orthofinder.sh`. The latter accepts the configured eleven-species
set without duplicating S_cotti, validates its frozen baseline checksum, and
refuses to reuse an existing run. Do not rerun either entry point merely to
synchronize code. These safeguards have synthetic regression coverage in
`scripts/test_entrypoint_safety.py` (66 tests across the full local suite).
