# *Schistocephalus cotti* Nuclear Annotation

This module records the *S. cotti* nuclear annotation independently from the
completed *S. solidus* workflow and from downstream comparative genomics.

## Server workspace

```text
/work/cyu/annotation/scotti-annotation
```

## Released baseline (20 November 2026)

The selected release is `/work/shared/Scotti_release_v1.0_preNCBI_20261118`:
9,234 genes, 11,464 protein isoforms and 9,234 longest representatives. See
[the authoritative release record](../release_v1.0.md).

Swiss-Prot, eggNOG and InterProScan have completed and their evidence is
integrated for retained IDs. The shared package includes the nuclear genome,
separate protein-coding and ncRNA GFF3 files, CDS and protein FASTAs,
mitochondrial working package, TE-review records and checksums. All three
shared GFF3 files passed validation. The complete baseline retains the
32 strong TE candidates and 47 unresolved cases; the 9,202-record exclusion
set is an optional sensitivity input only.

The 3,681 nuclear ncRNA loci remain candidates (344 tRNA + 3,337 Rfam).
Functional evidence is separate from curated product-name assignment. No
combined nuclear GFF3, release GTF or accessioned NCBI submission was generated.

## Historical search helpers

The scripts in this module retain the original BRAKER input paths and
`EXPECTED_PROTEINS=11471`, because the completed Swiss-Prot/eggNOG/InterPro
searches used all original isoforms. These are historical search parameters,
not the current release counts. Do not change them to 11,464 and overwrite
existing search prefixes merely to subset the retained proteins. Use the
published evidence tables and ID mappings.

Original BRAKER output: `results/braker3_etp_v1_retry2/` (9,241 genes,
11,471 isoforms). The filtered7 source and historical failed/retry logs remain
preserved. The later ncRNA/release scripts are maintained under
[`../post-mt/scripts/`](../post-mt/README.md).

## Layout

```text
scotti/annotation/
  config.example.sh
  scripts/
    common.sh
    01_validate_candidate.sh
    02_eggnog_mapper.sh
    03_interproscan.sh
    04_functional_status.sh
```

Copy and edit the server configuration:

```bash
cp scotti/annotation/config.example.sh scotti/annotation/config.sh
```

Run from the repository checkout on the server:

```bash
bash scotti/annotation/scripts/01_validate_candidate.sh
bash scotti/annotation/scripts/02_eggnog_mapper.sh
bash scotti/annotation/scripts/03_interproscan.sh
bash scotti/annotation/scripts/04_functional_status.sh
```

Search scripts reuse validated non-empty outputs. They do not overwrite the
BRAKER candidate or declare it to be a final release.

## Subsequent use

For comparative analysis use the published 9,234-record representative FASTA,
not all 11,464 isoforms and not the optional 9,202-record sensitivity set.
Existing frozen runs remain unchanged; verify input checksums before reusing
their results. New input choices belong in new run directories.

The background cestode species retain their public structural annotations and
receive standardized QC and functional annotation in `comparative-genomics/`.
Optional future model repair or NCBI submission does not invalidate the
documented v1.0 baseline.
