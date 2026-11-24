# Nuclear ncRNA, working snapshot and shared publication

This records completed server results supplied in the conversation through
20 November 2026. It does not launch new scans or reconstruct missing source data.

## Corrected nuclear tRNA selection

Input was the unmasked 412-sequence, 991,317,046-bp RagTag nuclear genome.
tRNAscan-SE 2.0.12 predicted 101,497 loci, many associated with repetitive
sequence. The original `--detail` run lacked the `--breakdown` HMM/secondary
structure columns required by EukHighConfidenceFilter. The resulting one-SeC
high-confidence result was a compatibility artefact and was superseded.

The corrected scan used `--detail --breakdown`. Its outputs are under
`results/ncrna_v1/01b_trnascan_breakdown/`; the corrected filter output is under
`03b_trna_highconfidence_breakdown/`. Official filter totals: 95,309 possible
pseudogenes, 5,840 secondary-cutoff failures, zero tertiary failures, 348 after
filtering, four isotype mismatches and **344 high-confidence loci**. The four
mismatches remain separate. The SeC locus is predicted, not functionally validated.

`04_trna_candidate_v1/` exports the 344 loci as flat precursor tRNA GFF3/BED
features with intron-column provenance and CDS-overlap audit. Mature sequences
and exon/intron hierarchies were not constructed.

## Rfam candidates

The completed `cmscan --cut_ga --rfam --nohmmonly --fmt 2 --clanin` scan has
31,773 raw hits in `02_rfam/Scotti_nuclear.rfam.tblout`. Rfam release/version
must be read from the actual database provenance; it is not inferred from the
solidus documentation. The scan was not repeated for candidate construction.

The project core policy retains included, nontruncated hits from explicitly
listed families with model-span coverage >=80%. The 80% threshold is a project
heuristic, not an Rfam requirement or proof of full transcript/function.
Original lower-scoring overlap mark `=` is excluded from the core; `*`, `^`
and `$` remain eligible. Residual same-strand candidate overlaps are resolved
by E-value, bit score and a deterministic ID tie-break. All raw hits/reasons
remain in the audit, including partial and non-core-family matches. Rfam tRNA
and tRNA-Sec matches are cross-checks and are not added to the 344-locus set.

`05_ncrna_candidate_v1/` contains **3,337 Rfam + 344 tRNA = 3,681 loci**.
GenomeTools passed. Zero Rfam candidates overlap same-strand protein CDS.
mir-71 is a precursor-homology candidate; mature miRNA coordinates were not assigned.

## Repeat-context audit

`06_ncrna_repeat_context_v1/` checked 1,353,789 RepeatMasker records. No candidate
or protein model was removed. Categories are independent interval unions;
Small RNA is not pooled into TE class and percentages need not add to 100%.

- TE-class overlap >=50%: 260 candidates; >=80%: 236.
- Small-RNA-class overlap >=80%: 657 candidates.
- U1: 970 candidates; zero TE-class or Small-RNA-class overlap reported.

These are review/context flags, not functional RNA counts or automatic
pseudogene calls. Zero overlap does not establish function. The cotti candidate
rules differ from the solidus conservative ncRNA construction, so the two
headline ncRNA counts are not directly comparable as biological copy numbers.

## Frozen source and selected release

Working snapshot: `/work/cyu/annotation/scotti-annotation/results/working_snapshots/Scotti_annotation_working_v1`.
It froze 63 files, checked gene/transcript/protein/CDS links and the exact
without32 subset, and validated protein/ncRNA GFF3. The genome was recorded by
path/hash at this stage. The 672 BUSCO group statuses matched exactly between
9,234 baseline representatives and the 9,202-record sensitivity set.

The user selected the **9,234-gene / 11,464-isoform / 9,234-representative
baseline** as scientific release v1.0. The 32 strong TE-associated candidate
genes and 47 unresolved cases are retained. The sensitivity set remains in
provenance and is not promoted to the release proteome.

Shared package: `/work/shared/Scotti_release_v1.0_preNCBI_20261118`.
It was published on **20 November 2026 at 19:03:34 UTC** with the same five
standard top-level component directories used by solidus, plus the preserved
working snapshot. Nuclear and mitochondrial genomes are included. Source/copy
checksums, three GFF3 validations and shared-group permissions passed.

The first publication stopped at an incorrect `gene == 36` check. The
mitochondrial GFF3 actually has **12 `gene` + 24 `ncRNA_gene`**, with 12 mRNA/CDS,
22 tRNA, two rRNA and 36 exon features. The corrected publisher checks the
combined locus types and supports checksum-verified `--resume`; it reused the
existing `.building` copy, then published the completed directory. Original
sources and existing comparative inputs were unchanged.

See [release files and validation](../release_v1.0.md). Builders under
`scripts/ncrna/` and `scripts/release/` preserve their original helper names.
The fixture tests use synthetic data; the publication fixture stubs external
GenomeTools/nuclear-size checks locally, while the real server publisher runs
those checks against the actual data.
