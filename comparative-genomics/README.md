# *Schistocephalus* Comparative Genomics

This directory contains the downstream multi-species comparative-genomics
workflow. It is intentionally separate from `../workflow-scripts/`, which
documents assembly and annotation of the *S. solidus* reference genome.

## Current species set

- `S_solidus`
- `L_intestinalis`
- `Spirometra_Aus1`
- `D_latus`
- `H_microstoma`
- `E_multilocularis`
- `E_granulosus`
- `T_solium`
- `T_multiceps`

The fragmented *S. erinaceieuropaei* PRJEB1202 annotation is retained only for
optional sensitivity analyses. The final analysis must be rerun after the
curated *S. cotti* annotation is available.

## Layout

```text
comparative-genomics/
  config.example.sh
  scripts/
    00_standardize_new_cyclophyllidea.sh
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

1. Standardize newly downloaded taxa and validate all inputs.
2. Calculate assembly/protein QC and run genome- and protein-mode BUSCO.
3. Run OrthoFinder and infer the species tree.
4. Run uniform InterProScan and eggNOG-mapper annotation.
5. Audit extreme gene families for TE and fragmented-gene artefacts.
6. Run DeepTMHMM2 and, when available, SignalP 6 for secretome analysis.
7. Add *S. cotti*, rerun the final orthology/tree analyses, and only then run
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

Combine the configured protein-mode BUSCO summaries with:

```bash
source comparative-genomics/config.sh
python comparative-genomics/scripts/02c_busco_summary.py \
  --busco-root "${QC}/busco_protein" \
  --species "${SPECIES[@]}" \
  --output "${QC}/all_species_protein_busco.tsv"
```

## Current completion boundary

The nine-species background dataset can be used for QC, functional annotation,
secretome analysis, exploratory OrthoFinder analyses, and preliminary gene
trees. Final OrthoFinder, CAFE5, focal-species gain/loss, dN/dS, branch-site,
TE-gene association, and microsynteny analyses must wait for the final
*S. cotti* annotation.

Low-completeness annotations, especially *L. intestinalis*, *D. latus*, and
*T. multiceps*, must not be used as sole evidence for gene loss or family
contraction.
