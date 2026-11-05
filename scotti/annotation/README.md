# *Schistocephalus cotti* Nuclear Annotation

This module records the *S. cotti* nuclear annotation independently from the
completed *S. solidus* workflow and from downstream comparative genomics.

## Server workspace

```text
/work/cyu/annotation/scotti-annotation
```

The current BRAKER3 ETP candidate is:

```text
results/braker3_etp_v1_retry2/
  braker.gtf
  braker.aa
  braker.codingseq
```

The candidate contains 11,471 protein predictions. It is not yet the final
primary proteome: structural and functional evidence must be integrated,
unsupported models reviewed, representative isoforms selected, identifiers
audited and a versioned release built.

## Current status (5 November 2026)

- Nuclear assembly and repeat masking are complete.
- RepeatMasker reports 69.93% interspersed repeats; this is annotation context,
  not a reason to discard genes by GFF keyword matching.
- BRAKER3 ETP produced the current 11,471-protein candidate set.
- Swiss-Prot search is complete.
- eggNOG-mapper is complete (`Total hits processed: 9225`).
- InterProScan was started successfully on all 11,471 sanitized proteins.
- Final evidence integration, model review, primary-isoform selection, release
  QC and publication of the release files are pending.

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

## Required completion order

1. Validate genome/GFF/GTF/CDS/protein relationships and sequence IDs.
2. Complete Swiss-Prot, eggNOG and InterProScan evidence.
3. Integrate functional evidence and review unsupported, very short, very long,
   internally stopped and TE-associated models.
4. Select one documented representative protein isoform per gene.
5. Run genome- and protein-mode BUSCO and OMArk where available.
6. Compare gene count, BUSCO, length distributions and annotation support with
   the *S. solidus* release.
7. Produce versioned genome, GFF3, CDS, all-isoform protein and primary-protein
   files with checksums.
8. Only after these checks, link the release into
   `02_standardized_data/` and add `S_cotti` to final comparative analyses.

The background cestode species are not reannotated here. They retain their
public structural annotations and receive standardized QC and functional
annotation in `comparative-genomics/`.
