# Comparative genomics status — 20 November 2026

This combines the user-reported completed analyses with the selected S. cotti
release. Updating local code/docs does not replace server inputs or launch jobs.

## Selected cotti release

`/work/shared/Scotti_release_v1.0_preNCBI_20261118` was published successfully
on 20 November 2026. The baseline is **9,234 genes/representatives and 11,464
protein isoforms**, retaining the 32 strong TE-associated candidate genes and
47 unresolved TE cases. The 9,202-record without32 set is optional sensitivity
provenance only. [Release details](../scotti/release_v1.0.md).

Use `proteins/Scotti_proteins.longest_isoform.v1.0.fa` for future baseline
input freezes. This publication did not alter the existing standardized
species links, frozen OrthoFinder inputs, orthogroups, trees or figures.
Compare existing input and release checksums before deciding to rerun.

## Eleven-species design and BUSCO

Focal species: S_solidus and S_cotti. Diphyllobothriidea backgrounds:
L_intestinalis, D_latus, Spirometra_Aus1 and S_proliferum. The five
Cyclophyllidea backgrounds are E_multilocularis, H_microstoma, E_granulosus,
T_solium and T_multiceps. The six-taxon core is not an instruction to remove
those backgrounds. S_erinaceieuropaei_PRJEB1202 remains excluded.

Reported protein-mode metazoa_odb12 results (same 672-group local dataset):

| Species | C | S | D | F | M |
|---|---:|---:|---:|---:|---:|
| S_solidus | 80.4 | 79.9 | 0.4 | 7.0 | 12.6 |
| S_cotti | 73.5 | 72.5 | 1.0 | 9.8 | 16.7 |
| L_intestinalis | 28.0 | 27.7 | 0.3 | 17.3 | 54.8 |
| D_latus | 20.4 | 20.4 | 0.0 | 39.0 | 40.6 |
| Spirometra_Aus1 | 59.5 | 54.2 | 5.4 | 15.6 | 24.9 |
| S_proliferum | 63.2 | 58.0 | 5.2 | 19.6 | 17.1 |
| E_multilocularis | 75.9 | 74.7 | 1.2 | 11.2 | 12.9 |
| H_microstoma | 72.3 | 71.3 | 1.0 | 11.8 | 15.9 |
| E_granulosus | 76.3 | 74.3 | 2.1 | 9.2 | 14.4 |
| T_solium | 65.8 | 65.6 | 0.1 | 15.0 | 19.2 |
| T_multiceps | 57.1 | 56.4 | 0.7 | 17.3 | 25.6 |

Preserve odb10 for historical/model investigations; percentages from different
lineages are not changes in quality. The cotti baseline and without32 protein
sets have zero changes in all 672 BUSCO group statuses. Low-completeness
annotations, especially Ligula/D. latus, confound absence and contraction claims.

## Completed S. proliferum work

SprJ_4.0, GCA_053814205.1 / PRJDB8966: 29,231 genes/mRNAs, 28,993 translated
proteins, 238 no-CDS models and no multi-mRNA genes. All translated proteins
map to CDS IDs; no-CDS models are audited, not automatically pseudogenes.
InterProScan completed: 26,006 proteins with member-database matches and
18,640 with integrated InterPro entries (64.29%). eggNOG reported 21,019
processed hits. Do not infer a new run or secretome completion from these counts.
Recovery/mapping of the paper's reported TE-associated ID list remains separate.
TE exclusion is an explicit analysis choice, not a mandatory OrthoFinder step.

## Completed OrthoFinder baselines

Two separate 11-species runs used the 9,234-record filtered7 cotti baseline:

- Original auto-rooted: `04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/OrthoFinder/Results_Nov09/`.
  Final membership/count tables agree at **16,583 orthogroups / 137,611 assigned
  proteins**. Earlier statistics reported 16,669 groups / 137,591 proteins;
  those summary values must not replace the final membership counts.
- Full rerun with the supplied Cyclophyllidea-rooted tree:
  `04_orthofinder/runs/exploratory11_cyclophyllidea_rooted_full_v1/Results_Nov09/`.
  Final membership/count tables agree at **13,937 groups / 138,538 assigned
  proteins**, with no duplicate protein memberships. The run summary reported
  13,975 / 138,533; membership-based figures use the final tables.

These runs were not uniformly TE-filtered. Publication of the same cotti
baseline does not automatically change their inputs or remove these analysis
limitations. A whole-proteome without32 run is optional, not the selected
release or a prerequisite to acknowledging the published baseline.

Local historical records mark `exploratory11_scotti_without32_v1` invalid:
cotti overwrote the staged T_multiceps file, leaving ten inputs. Script 16
fixes staging and targets `_v2`; successful completion of that corrected
run is not established by the outputs supplied here. Keep the invalid run
out of phylogenetic/family interpretations.

## Independent phylogeny and figure boundary

`08_phylogeny/singlecopy374_v1` contains **373**, not 374, strict single-copy
loci selected from the final membership table and matched to the corresponding
alignments. The cached 374-group list/single-copy sequences are inconsistent
with those final group IDs and must not be mixed with them. The supermatrix
has 190,437 amino-acid sites and 11 species.

IQ-TREE support values are SH-aLRT/UFBoot percentages, not OrthoFinder node
support. The short deep edge has gCF 38.61% and sCF 31.89% despite 100/100
support. Gene-tree resolutions: 144/373 Sparganum/Spirometra basal, 68/373
Schistocephalus basal, 129/373 D. latus/Ligula basal, 32 other discordances.
The cotti/solidus pair is stable (gCF 86.06%, sCF 92.20%). Rooting between
the sampled orders is an explicit outgroup choice.

Do not transfer internal duplication counts from the OrthoFinder topology
onto the different IQ-TREE topology. Terminal duplication summaries may be
shown separately with their source stated. The historical plot script 11
reproduces the earlier auto-rooted figure; changing release labels does not
recompute its embedded values or fix its historical topology.

## Configuration and next use

The default SPECIES array remains the ten backgrounds to preserve existing
functional/secretome/MEROPS jobs. S_cotti is released, not pending; script 03
can explicitly include the shared release with `INCLUDE_SCOTTI_RELEASE=1`.
Historical candidate mode remains available for reproducing earlier input
freezes. The two modes are mutually exclusive. No standardized links or
server configurations were changed by this local update.

Uniform family-count interpretation still requires auditing extreme/TE-related
families and annotation completeness. CAFE5 and molecular-evolution analyses
need analysis-specific validated inputs. Existing secretome/MEROPS scripts
and computational parameters are preserved.
