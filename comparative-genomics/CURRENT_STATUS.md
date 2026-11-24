# Comparative genomics analysis status

Updated 24 November 2026, America/Chicago. This is a record of supplied server
results and locally inspected files, not a new server audit. Code available,
output present, and analysis validated are different states. The local update
does not submit jobs, replace releases, or modify previous result directories.

## Analysis overview

| Analysis | Recorded state | Next verification |
|---|---|---|
| Selected nuclear/protein releases | S. solidus 0831 and S. cotti 1118 selected; cotti 9,234-protein release/input SHA256 agrees | Preserve the baseline; accessioning is separate |
| Eleven-species OrthoFinder baseline | Completed; Figure 1 uses the 9,234-protein cotti input | Reconcile overall Statistics versus final membership totals before reporting global OG/assigned totals |
| without32 sensitivity | Valid eleven-species `_v2` completed; `_v1` had only ten species and is invalid | Do not substitute sensitivity statistics for baseline Figure 1 |
| Baseline phylogeny | 373 strictly single-copy loci, 190,437 aligned amino-acid sites | Retain deep-node discordance caveat |
| Genome and protein BUSCO | All 22 selected odb12 assessments reported complete | Preserve input identities and matching dataset |
| Table 1 v4 | Collected; D. latus exon hierarchy overcount diagnosed | Corrected collector tested locally; v5 server collection not yet confirmed |
| Old S. solidus assemblies | Two assemblies downloaded, assembly statistics and genome BUSCO reported | Keep these as within-species assembly comparisons, not two added species |
| Figure 1 | Revised real-data PDF/SVG/PNG produced locally | Use 16/20 baseline, larger labels, dashed extensions and marked branch breaks |
| Raw focal species-specific OGs | 16 cotti / 20 solidus extracted; sister-protein/genome searches completed | Raw clustering is not proof of novel family origin |
| Prioritized species-specific candidates | RNA support reported for 15 proteins; cotti two-protein CDS repeat overlap is zero | Solidus candidate repeat overlap, mapping specificity and model validity still need review |
| Functional annotation | Focal release evidence integrated; S. proliferum InterPro/eggNOG completion reported | Uniform eleven-species coverage and ID/version matching need current output audit |
| Eight-species CAFE5 screening | Gamma-k3 fits v2/v3/v4 reported with near-identical likelihoods | Treat calibration, annotation and family definitions as screening limitations |
| Figure 2 | User supplied completed all-change tree/pie figure | Counts are not significance-filtered; current server figure hashes not retrieved |
| Four CAFE priority families | Initial annotations/copy counts available; script 31 format fix available | Successful post-fix per-protein review and RNA/repeat/gene-tree validation not yet reported |
| SignalP / DeepTMHMM2 / secretome / MEROPS | Workflow scripts available; eleven-species completion not established in this record | Inventory existing outputs before any rerun |
| dN/dS / branch-site / family trees / microsynteny | No completed analysis documented here | Define validated inputs; reference-guided scaffold order is not independent synteny evidence |

## Selected inputs and version boundaries

Server project: `/work/cyu/Schistocephalus_comparative_genomics`.

- S. solidus release: `/work/shared/Ssolidus_release_v1.0_preNCBI_20260831`.
- S. cotti release: `/work/shared/Scotti_release_v1.0_preNCBI_20261118`.
- Baseline: `/work/cyu/Schistocephalus_comparative_genomics/04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/OrthoFinder/Results_Nov09`.
- Sensitivity: `/work/cyu/Schistocephalus_comparative_genomics/04_orthofinder/runs/exploratory11_scotti_without32_v2/input/OrthoFinder/Results_Nov12`.

The selected cotti release, frozen working snapshot, and baseline OrthoFinder
input contain the same 9,234 representatives, with SHA256
`07cf93ad5b6edf2c27e5d42745fe425cd9b6029af1810aaf63ab8ee2b0d1d371`.
The 32 reviewed strong-TE candidates remain in this release; 9,202 proteins
belong to the separate without32 sensitivity input. There is no uniform
eleven-species TE-protein census. The 47 unresolved cotti TE-review cases were
not automatically removed.

The selected cotti assembly is RagTag-only after removing four whole
mitochondrial candidate contigs, not either SAMBA candidate. It contains
412 sequences and 991,317,046 bp. Its scaffold N50 reflects scaffolding
against solidus. The original complete RagTag command/AGP and whether
`ragtag correct` was used are not established by the retained local records.
Pre-scaffolding nuclear contig N50 and independent chromosome-scale
rearrangement evidence remain separate checks. The retained assembly record
is `scotti/post-mt/execution_record.md` in the main local workspace; that
upstream module is not included in the comparative-only synchronization pack.

## BUSCO and Table 1

The selected dataset is metazoa_odb12, creation date 2026-05-22, n=672.
Genome and representative-protein modes are separate assessments, not two
lineage versions. odb10 remains historical and must not be merged numerically.

| Species | Representatives | Genome C (%) | Protein C (%) |
|---|---:|---:|---:|
| S_solidus | 9,796 | 71.6 | 80.4 |
| S_cotti | 9,234 | 71.9 | 73.5 |
| L_intestinalis | 27,256 | 67.9 | 28.0 |
| D_latus | 19,930 | 26.0 | 20.4 |
| Spirometra_Aus1 | 11,064 | 73.7 | 59.5 |
| S_proliferum | 28,993 | 72.3 | 63.2 |
| E_multilocularis | 10,656 | 69.9 | 75.9 |
| H_microstoma | 10,139 | 67.1 | 72.3 |
| E_granulosus | 9,985 | 70.1 | 76.3 |
| T_solium | 11,854 | 69.6 | 65.8 |
| T_multiceps | 12,502 | 71.4 | 57.1 |

Total representatives: 161,409. Complete batch record:
`/work/cyu/Schistocephalus_comparative_genomics/03_qc/final_genome_busco_odb12_v1`.
The preserved table is
`/work/cyu/Schistocephalus_comparative_genomics/03_qc/release_table1_11species_v4/Table1_11species.tsv`.

D. latus has 111,669 raw exon rows: 55,608 linked to genes, 55,835 to mRNAs,
and 226 to pseudogenes. The transcript-only audit reported **55,835 total
associations, 55,608 coding associations, 227 pseudogene-transcript
associations, 2.79 exons per coding transcript, and 227.9 bp mean coding exon
length**. The extra hierarchy rows are not grounds to delete source GFF data.
Collector 30 now uses `transcript_exon_associations_v2`, with per-species
`exon_transcript_audit.tsv`. The new server destination is
`/work/cyu/Schistocephalus_comparative_genomics/03_qc/release_table1_11species_v5_exonfix`;
no completion output from that collection has been supplied yet.

E. multilocularis has 9,192 explicitly linked coding gene IDs, 10,656 protein
IDs, and 1,464 proteins whose CDS lacks a linked gene feature. Its GFF has no
exon features: do not invent exon metrics or equate all protein records with
explicit gene features. V4's 24 missing cells are 11 TE-protein counts, nine
background RepeatMasker percentages and four E. multilocularis exon metrics.
V5 missing-cell counts must be inspected rather than assumed unchanged.

Solidus repeat-masked bp are 580,875,018 / 846,268,314 (68.63958%); both
numerator and denominator exclude the 6,055-bp fish scaffold_61. Retained
sequences were verified. Cotti has 715,044,557 / 991,317,046 (72.13076%);
the supplied audit established coordinates/lengths, not sequence identity
against the original RepeatMasker input. That last identity check remains open.
Neither result is a TE-protein count.

Old solidus comparisons: GCA_900618435.1, 539,424,504 bp, 56,778 sequences,
N50 31,582 bp, genome BUSCO C61.9%; GCA_017591395.1, 625,207,408 bp,
15,357 sequences, N50 121,651 bp, C67.4%. The matched workflow used
BUSCO/Miniprot and the same odb12 lineage. Internal-stop warnings were
reported; these scores do not establish absence of gene-model issues.

## Figure 1 and phylogeny

The baseline phylogeny uses 373 loci / 190,437 sites, whereas without32 v2
uses 381 / 197,179. Both have eleven species. Keep them labelled separately.
The baseline's short deep split has gCF 38.61% and sCF 31.89% despite
SH-aLRT/UFBoot 100/100; the cotti/solidus split has gCF 86.06% and sCF 92.20%.
Rooting between the sampled orders is explicit. Terminal duplications are
separate OrthoFinder summaries, not CAFE family changes or counts reconciled
onto every branch of the displayed IQ-TREE tree.

| Species | Terminal duplications | Assigned proteins | Species-specific OGs |
|---|---:|---:|---:|
| S_solidus | 628 | 9,540 | 20 |
| S_cotti | 692 | 8,947 | 16 |
| L_intestinalis | 11,145 | 19,974 | 1,453 |
| D_latus | 4,520 | 16,232 | 141 |
| Spirometra_Aus1 | 1,836 | 10,492 | 55 |
| S_proliferum | 9,138 | 21,161 | 1,001 |
| E_multilocularis | 865 | 9,956 | 54 |
| H_microstoma | 1,614 | 9,355 | 206 |
| E_granulosus | 756 | 9,470 | 35 |
| T_solium | 1,560 | 10,757 | 35 |
| T_multiceps | 3,813 | 11,707 | 221 |

The final local Figure 1 is
`output/pdf/Figure1_baseline9234_largefont_dashed_v3/orthology_summary_clean.pdf`
under the repository root, with matching SVG, PNG, source data and branch-break
audit. Only Spirometra is italicized in "Spirometra sp. Aus1"; Sparganum
proliferum remains unabbreviated. Dashed guides are not branch lengths, and
double-slash marks identify display-compressed long branches. Old 10/13
species-specific counts and 705/643 duplication labels are superseded.

Global baseline membership summaries reported 16,583 OGs / 137,611 assigned
proteins / 1,460 all-species OGs; the saved local Statistics_Overall reports
16,669 / 137,591 / 1,445. Do not combine these totals or claim that the
discrepancy has been resolved. The focal Figure 1 values and 373 strict loci
are the separately checked values above. A separate Cyclophyllidea-rooted
full OrthoFinder rerun is another result, not a replacement by relabelling.

## Species-specific candidate evidence

The raw baseline contains 16 cotti-only OGs / 74 proteins and 20 solidus-only
OGs / 97 proteins. Sister-proteome BLASTP plus sister-genome TBLASTN left
seven cotti and eighteen solidus proteins unresolved at the chosen coverage
screen. Lack of a retained hit is not proof of biological absence.

RNA/repeat priority membership is **two cotti plus thirteen solidus proteins**,
not the four-family CAFE priority set below. Cotti OG0015711 has g8804.t1
(138 aa, strong RNA coverage) and g8806.t1 (226 aa, partial/moderate coverage);
both have zero measured CDS repeat overlap. Their release coordinates are on
the same strand of scaffold_8_RagTag, within 94.9 kb. Whether they lie within
one original pre-RagTag contig has not been verified; do not treat that
scaffold spacing as independent tandem-duplication evidence.

Solidus OG0016342 (two proteins), OG0016337 (two), and OG0009067 (six) have
reported strong RNA coverage. OG0012892 (three) has weak/partial support.
Solidus candidate repeat overlap and uniqueness of RNA mapping are still
separate validation tasks. Transcription does not establish translation or
species specificity. OG0002484 has Fibronectin type III-containing proteins
and homologs in cotti; it is not a demonstrated newly originated family.

## Functional annotation and other downstream analyses

The focal releases include integrated functional evidence. For S. proliferum,
InterProScan completion was reported with 26,006 any-member matches (89.7%),
18,640 integrated InterPro matches (64.3%) and 16,144 Pfam matches (55.7%)
among 28,993 proteins. eggNOG completed with 21,019 non-comment annotation
rows; that row count is not a cross-species uniform-coverage assessment.
The documented 29,231 gene/mRNA models include 238 untranslated models.

All-eleven functional coverage, SignalP, DeepTMHMM2, integrated secretome,
MEROPS, dN/dS, branch-site selection and family-tree completion need actual
output evidence. Scripts or old subset results do not establish completion
for the current eleven-species inputs. No new searches were launched here.

## CAFE5 screening and Figure 2

The analysed eight species are S_solidus, S_cotti, Spirometra_Aus1,
E_multilocularis, H_microstoma, E_granulosus, T_solium and T_multiceps.
Ligula/D. latus were excluded for low protein completeness; S. proliferum
was omitted from this screen pending analysis-specific TE-family review.
The eight-species screen is not a completed, uniformly filtered eleven- or
nine-species final analysis. The planned species-set variants have not been
reported complete here.

The baseline cotti TE32 list mapped 27 proteins to 18 excluded OGs; five
were unassigned and already outside the family matrix. Matrix preparation
retained 10,480 families, separated two large families, and excluded 6,101.
CAFE root-presence filtering retained 6,628 and removed 3,852.
Time-tree calibration fixed the T. solium/T. multiceps split at an
analyst-selected secondary age of 16.5 Ma, producing a reported root height
of 54.19505613 Ma. Its source topology was the without32 v2 IQ-TREE tree;
the family counts came from the 9,234 baseline. This mixed provenance is
explicit, not an assertion that both analyses used the identical input set.

Run root:
`/work/cyu/Schistocephalus_comparative_genomics/07_gene_family/cafe5/runs`.

| Gamma-k3 fit | -lnL | Lambda | Alpha | Epsilon |
|---|---:|---:|---:|---:|
| screening8_16.5Ma_k3_v2 | 39975.4 | 0.0047048280954719 | 0.375730 | 0.0219204 |
| screening8_16.5Ma_k3_v3 | 39975.4 | 0.0047044326550371 | 0.375788 | 0.0219204 |
| screening8_16.5Ma_k3_v4 | 39975.4 | 0.0047043757458531 | 0.375674 | 0.0219204 |

These fits show similar reported numerical optima, not proof of global
optimality or biological robustness. The initial v2 FAILED flag referred to
the validator expecting Gamma_change.txt instead of the actual .tab file;
the subsequent output check passed. V1 had a separate argument-parsing error
and is not the completed Gamma analysis.

| Terminal lineage | All expanded | All contracted | Nominal family-P<0.01 expanded | Nominal family-P<0.01 contracted |
|---|---:|---:|---:|---:|
| S_cotti | 265 | 611 | 14 | 39 |
| S_solidus | 178 | 230 | 24 | 18 |

All-change counts for the other terminal lineages are E_granulosus +241/-393,
E_multilocularis +202/-124, T_multiceps +194/-1360, T_solium +718/-251,
H_microstoma +404/-674, and Spirometra_Aus1 +644/-935. The branch leading to
the cotti/solidus ancestor has +75/-226. These count changing families, not
total gene copies gained or lost.

Figure 2 shows **all changes among 6,628 families**, without a significance
filter. Internal pies refer to ancestral incoming branches; tip pies refer
to terminal branches. Blue is expansion, vermillion contraction, grey unchanged.
Do not label these counts "significantly expanded/contracted".

The nominal counts were reported in supplied summary output. They are
family-wide exploratory screens, not branch-specific tests. Script 29 schema
v2 keeps Gamma branch values as reported posterior probabilities, never as
P-values or candidate gates. Old familyFDR_branchP / Branch_pvalue_raw results
are retired. A successful server execution of the latest schema-v2 summary
has not yet been supplied; use a new output directory when verifying it.

Reported family-P diagnostics: 106 zeros, smallest positive value 0.001, and
316 raw P<0.05 families. Replacing values below 0.001 by that floor for a
separate BH sensitivity calculation produced no q<0.05 expanding or contracting
families for either focal species in the supplied summary. This is a
resolution-sensitivity assumption, not an estimate of the true P-values or
evidence that no biological changes occurred. Do not promote q=0 from reported
P=0 into a definitive FDR claim.

### Four CAFE priority families

| OG | Initial domain description | Cotti copies | Solidus copies |
|---|---|---:|---:|
| OG0000875 | Tctex-1-like | 9 | 1 |
| OG0000880 | CHORD/CS | 8 | 2 |
| OG0000999 | Ribonucleotide reductase small subunit | 7 | 1 |
| OG0001937 | CRAL-TRIO | 2 | 7 |

These are priority candidates, not confirmed branch-significant expansions.
Script 31 targets 26 cotti and 11 solidus proteins. Its earlier g4713.t1 MD5
stop was traced to a single terminal '*': release/OrthoFinder FASTA stores
81 characters, while the InterPro search used the same 80-aa sequence without
the marker. Comparison-only normalization is implemented and tested; release
and analysis FASTAs remain unchanged. The complete post-fix review output
has not been reported, nor have the CAFE-family RNA/repeat and gene-tree
checks. Do not reuse the unrelated fifteen-protein candidate validation as
evidence that these thirty-seven proteins have been reviewed.

## Synchronization and next checks

### Manuscript entry-point corrections (24 November 2026)

The selected releases, baseline Figure 1 counts/alignment, and existing
eight-species CAFE5 fits are unchanged. The old 10-species summary default is
corrected to eleven species / odb12. Staging separates 11,464 all-isoform and
9,234 representative cotti proteins and refuses conflicting destinations.
InterPro runs use derived terminal-stop-normalized input and verifiable
completion records; historical results without these records are preserved,
not automatically adopted or rerun. The recovered TagSeq script now handles
empty exclusion lists safely, checks GNU awk, and refuses existing final GTFs.
Historical biological outputs were not regenerated or retrospectively certified.

The user-selected main local VS Code tree is
`/Users/changyueyu/Desktop/Worm-genome-assembly`.
The copies under
`/Users/changyueyu/Documents/Codex/2026-08-25/ta/Worm-genome-assembly` and
`/Users/changyueyu/Documents/Codex/2026-10-01/hi-kirsty-i-d-like-to/outputs/Worm-genome-assembly`
are secondary historical/synchronization copies, not automatically current mirrors.
The Desktop solidus workflow now includes the recovered TagSeq script 21 and
Rfam filtering additions to scripts 12–13. No biological analysis was rerun.
The cotti mitochondrial v2 builder was subsequently copied by the user from the
server release package to `scotti/mitochondria/scripts/09_build_curated_candidate_v2.py`.
Local syntax validation passed; the script requires Biopython and the original
server MITOS2 GFF/unrotated FASTA. No curation was rerun and a remote checksum
comparison remains unrecorded. Historical local script 08 is not its replacement.
This edit does not establish synchronization with the server or GitHub.
`scripts/36_audit_comparative_outputs_readonly.py` inventories expected server
paths and can compare code hashes against `data/analysis_sync_manifest.json`.
Presence checks are explicitly not full completion validation.

1. Confirm the corrected v5 exon collection and inspect every species' audit.
2. Verify the server's corrected CAFE summary schema and resume the four-family
   review into a new output directory; reuse existing CAFE fits/searches.
3. Audit the existing eleven-species functional/secretome/MEROPS outputs before
   planning missing work. Finish source/version, repeat-input identity and
   orthogroup-total discrepancies before final publication claims.

### GitHub synchronization preparation (2026-11-24)

Compared this source set against the live GitHub main branch before staging.
Preserved the upstream BUSCO lineage-specific output directories and restored
OrthoFinder's fresh-run input copies, checksum recording and refusal to reuse
existing runs. The OrthoFinder entry point accepts S_cotti already in the
configured eleven-species set and verifies the selected baseline SHA256.
All 66 synthetic/local regression tests passed; no biological analysis was
rerun and no server inputs, results or completion records were changed.
