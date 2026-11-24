# S. cotti scientific release v1.0 (pre-NCBI)

Publication was confirmed by user-supplied stickleback output on **20 November 2026 at 19:03:34 UTC**. The directory's `20261118` suffix is the selected version label, not the completed copy date. This local update does not submit server jobs or modify published data.

## Authoritative shared package

```text
/work/shared/Scotti_release_v1.0_preNCBI_20261118
```

| Component | Selected release |
|---|---:|
| Nuclear sequences | 412 |
| Nuclear bases | 991,317,046 |
| Protein-coding genes | 9,234 |
| Protein isoforms / CDS sequences | 11,464 |
| Longest representative proteins | 9,234 |
| Nuclear ncRNA candidate loci | 3,681 (344 tRNA + 3,337 Rfam) |
| Mitochondrial bases | 15,437 |
| Mitochondrial genes | 36 (12 PCG + 22 tRNA + 2 rRNA) |

The user selected the complete filtered7 baseline. All 32 strong TE-associated candidate genes and all 47 unresolved TE cases remain in the release. The 9,202-record without32 FASTA is a sensitivity input retained only in provenance; it is not the selected release proteome.

## Files relative to the shared package

| Path | Purpose |
|---|---|
| `assembly/Scotti_nuclear.fa` | Selected RagTag nuclear genome, four whole mt candidates removed |
| `assembly/Scotti_mitochondrial.fa` | Separate mitochondrial genome |
| `assembly/Scotti_mitochondrial.candidate.v1.0.gff3` | Curated mitochondrial candidate annotation |
| `assembly/mitochondria_working_v1/` | Complete preserved mitochondrial working package, including the authoritative server v2 builder |
| `annotation/Scotti_nuclear.protein_coding.v1.0.gff3` | Released filtered7 protein-coding annotation |
| `annotation/Scotti_nuclear.CDS.all_isoforms.v1.0.fa` | All 11,464 CDS sequences |
| `annotation/Scotti_nuclear.ncRNA.candidate.v1.0.gff3` | Separate 3,681-locus ncRNA candidate component |
| `proteins/Scotti_proteins.all_isoforms.v1.0.fa` | All 11,464 protein isoforms |
| `proteins/Scotti_proteins.longest_isoform.v1.0.fa` | Selected 9,234 representatives for protein comparisons |
| `functional_annotation/` | ID-linked Swiss-Prot, eggNOG and InterPro evidence |
| `working_snapshot/TE_sensitivity/` | Review decisions and the 9,202-protein sensitivity set |
| `working_snapshot/qc/BUSCO_odb12/` | Baseline/without32 BUSCO tables and summaries |
| `quality_control/` | Source/copy checksum logs and three GFF3 validation logs |
| `file_map.tsv`, `summary.json`, `STATUS.txt`, `SHA256SUMS` | Release mapping, status, counts and checksums |

Protein-coding and ncRNA GFF3 files remain separate. No combined nuclear GFF3 or release GTF was constructed. Functional tables are evidence, not curated product-name assignments. The mitochondrial and ncRNA components retain their documented candidate status; release publication does not establish every ncRNA's functional or pseudogene status.

## Confirmed checks

The source and copied packages passed checksums. Nuclear protein-coding, nuclear ncRNA and mitochondrial GFF3 files all passed GenomeTools validation. Shared-group permissions were applied; original `/work/cyu` data were preserved. The package is self-contained for the included files; standard-layout aliases link only to newly copied files within the shared package.

Nuclear FASTA SHA256:

```text
e51d84e4d3be27821474b0550f394376533f111fe14eb219ce9c054c786bc9a6
```

BUSCO metazoa_odb12 (672 groups): `C:73.5%[S:72.5%,D:1.0%],F:9.8%,M:16.7%`. This is 494 complete groups: 487 single-copy plus 7 duplicated, with 66 fragmented and 112 missing. The baseline and without32 sensitivity sets have **zero per-group status changes**, not merely matching aggregate percentages.

The ncRNA repeat-context audit flagged 236 candidates with TE-class overlap >=80%; these were retained for review. All 970 U1 candidates had zero reported TE-class and Small-RNA-class overlap. Absence of such overlap does not prove function. The 657 high Small-RNA-class overlaps are not automatically TE or pseudogene calls.

## Reproducibility and use

- [ncRNA, snapshot and publication record](post-mt/ncrna_and_release.md)
- [Historical filtered7 source package](post-mt/filtered7_working_candidate.md)
- [Comparative analysis status](../comparative-genomics/CURRENT_STATUS.md)

Local publication code is under `post-mt/scripts/release/`. The corrected publisher counts 12 `gene` plus 24 `ncRNA_gene` loci and supports explicit `--resume` of a checksum-matching `.building` directory. Publication completed; do not rerun the publisher against the existing release.

Existing comparative inputs and results were not replaced. Compare a frozen OrthoFinder input with the released representative FASTA by checksum before deciding whether a new run is necessary. Renaming/copying identical proteins into a release does not itself require repeating orthology, functional searches or phylogeny.

This is a scientific pre-NCBI release. Submission identifiers and NCBI submission/validation are a separate stage.
