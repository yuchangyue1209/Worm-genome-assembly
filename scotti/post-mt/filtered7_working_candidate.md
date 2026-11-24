# Historical filtered7 source of the selected v1.0 release

Original candidate-stage record: 5 November 2026. This exact baseline was later
selected for the scientific v1.0 release, published on 20 November 2026.
See [the release record](../release_v1.0.md) for authoritative shared paths.
The historical source directory and its status labels are preserved.

## Paths and identifiers

Server package: `/work/cyu/annotation/scotti-annotation/results/filtered_candidates/Scotti_filtered7_v1`.

|File|Purpose|
|---|---|
|Scotti_nuclear.filtered7.candidate.gff3|Protein-coding candidate annotation|
|Scotti_nuclear.filtered7.cds.fa|11,464 coding sequences|
|Scotti_nuclear.filtered7.proteins.fa|11,464 protein isoforms|
|Scotti_nuclear.filtered7.longest.proteins.fa|9,234 representative proteins|
|excluded_models.tsv; excluded.proteins.fa|Audit trail for excluded predictions|
|longest_isoform_selection.tsv|Representative transcript-to-gene mapping|
|summary.json|Candidate counts/status|
|functional_evidence_v1/|Filtered source TSVs, all/primary evidence matrices, coverage and checksums|

Genome unchanged: `/work/cyu/assembly/worm/scotti/final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa`, 412 sequences, 991,317,046 bp. BRAKER used its softmasked counterpart under repeats/02_repeatmasker. The mt genome remains separate.

## Filter decision

Excluded current models: g1218 (10 aa), g196 (5 aa), g3336 (8 aa), g5529 (10 aa), g5573 (4 aa), g7051 (4 aa), g7437 (2 aa). All lack Swiss-Prot, eggNOG and member-database/InterPro support.
These are exclusions of extremely short protein-coding predictions, not claims that no real gene exists at their loci.
- g196 lies 88 bp from same-strand g195.
- g3336 lies 61 bp from same-strand g3337 and has a 22-bp annotated CDS.
- g5529 lies 210 bp from same-strand g5528 and has substantial regional RNA signal.
- g7051 lies inside the span of same-strand g7050.
These loci remain available for repair review; no models were merged automatically. The 37 other short-gene candidates and two short alternative isoforms in longer genes were retained.

Before filtering: 9,241 genes / 11,471 isoforms. After: 9,234 / 11,464. Original BRAKER files are unchanged.

## Validation

GenomeTools reported `input is valid GFF3` after an initial header-only repair. Script 27 now ensures ##gff-version 3 is the first line. Missing ##sequence-region warnings were nonfatal and did not rewrite the file. Full FASTA coordinate and translation validation remain distinct release checks.
BUSCO metazoa_odb12: C73.5% [S72.5%, D1.0%], F9.8%, M16.7%, n672. Counts: C494, S487, D7, F66, M112. These counts are identical before and after filtering. This establishes no loss by this assay, not complete biological validation.
The initial SHA256SUMS predates the header repair; use a subsequently verified manifest. SHA256SUMS.validated creation was requested, but its successful check output has not been reported. Script 30 can create an additional dated top-level manifest. Functional evidence has its own manifest.

## Completed evidence integration

InterProScan 5.78-109.0 completed with `100% done` from approximately 16:20 to 18:11 on Nov 5, 2026. The original search used all 11,471 proteins. Filtered TSVs were selected by retained IDs without rerunning searches. The original large InterProScan GFF3 was not rewritten.

|Evidence|All 11,464 isoforms|9,234 representatives|
|---|---:|---:|
|Swiss-Prot hits|8,206|6,545|
|eggNOG annotation rows|9,225|7,351|
|InterPro member-database matches|10,594|8,490|
|Integrated InterPro entries|9,391|7,499|
|Any of the three sources|10,721|8,592|

642 representative proteins lack hits across these resources and are retained. Hits/annotation rows do not establish validated functions. Member signatures may be nonspecific. No automatic product names were assigned and the GFF3 remains structural.

## Release and comparison boundary

The selected scientific release retains all 9,234 representative proteins and
11,464 isoforms from this source. TE case review created a separate
9,202-record sensitivity set, but the user did not select it as the release.
The later snapshot verified matching baseline/without32 statuses for every
BUSCO group. Source/copy checksums and three shared GFF3 checks passed at
publication; use the shared package's top-level SHA256SUMS rather than an
older candidate checksum manifest that predates the header repair.

Nuclear ncRNA annotation is now available as a separate 3,681-locus candidate
component. Functional evidence is integrated; curated product-name assignment
and accessioned NCBI submission remain separate from this scientific release.
Existing comparative input files and results were not replaced. Cross-species
annotation completeness, TE-associated families and important candidate loci
still require consideration in interpreting gene-family changes.
