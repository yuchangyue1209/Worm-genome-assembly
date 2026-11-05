# *Schistocephalus* Comparative Genomics

This directory contains the downstream multi-species comparative-genomics
workflow. It is intentionally separate from `../workflow-scripts/`, which
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

These ten species form the current preliminary dataset. *S. cotti* is pending
and will bring the final dataset to eleven species. The fragmented
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
7. Recover the published *S. proliferum* TE-associated locus list, build a
   documented non-TE protein set and use it for final gene-family inference.
8. Add *S. cotti*, rerun the final orthology/tree analyses, and only then run
   the final CAFE5 and molecular-evolution analyses.

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

The ten-species preliminary dataset can be used for QC, functional annotation,
secretome analysis, exploratory OrthoFinder analyses, and preliminary gene
trees. Final OrthoFinder, CAFE5, focal-species gain/loss, dN/dS, branch-site,
TE-gene association, and microsynteny analyses must wait for the final
*S. cotti* annotation.

The public *S. proliferum* SprJ_4.0 annotation contains 29,231 gene/mRNA
models and 28,993 translated proteins. Each gene has at most one mRNA, so no
longest-isoform selection is required. The 238 untranslated models are retained
in an audit table but are absent from protein-based analyses. The paper reports
6,316 TE-associated loci; these must be recovered from the official
supplementary data rather than inferred by keyword matching the GFF3.

Low-completeness annotations, especially *L. intestinalis*, *D. latus*, and
*T. multiceps*, must not be used as sole evidence for gene loss or family
contraction.
