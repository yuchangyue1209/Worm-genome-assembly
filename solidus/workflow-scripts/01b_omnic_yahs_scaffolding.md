# Omni-C scaffolding with YaHS

This stage follows ONT assembly, Racon/Medaka polishing, and purge_dups,
and precedes mitochondrial-sequence separation and downstream annotation.

## Documented procedure

Omni-C reads were aligned to the purged assembly using BWA-MEM2. The
alignments were used for scaffolding with YaHS using the q10 and r500k
settings, with contig error correction disabled.

The source record does not include the original BWA-MEM2 or YaHS command
lines, software versions, alignment-filtering details, or Omni-C library
preparation information. The q10/r500k labels are preserved as recorded;
this document does not reconstruct unverified command-line flags.

## Assembly curation

Mitochondrial sequence was maintained separately from the nuclear assembly.
Twenty-six short, high-depth rDNA-associated contigs (82,735 bp in total)
were archived separately as unplaced repetitive sequences and were not
classified as exogenous contamination.

The subsequent FCS-GX audit identified the 6,055-bp fish-derived
scaffold_61 for exclusion. The final nuclear release comprises 166 scaffolds
and 846,268,314 bp, with a scaffold N50 of 99,718,931 bp.

See `03_document_mt_removal.sh` for mitochondrial separation and
`09_fcs_gx_contamination_audit.sh` for the final contaminant audit.

## Manuscript Methods text

Following long-read assembly, consensus polishing, and haplotig purging,
Omni-C reads were aligned to the assembly using BWA-MEM2. Scaffolding was
performed using YaHS with the documented q10 and r500k settings and contig
error correction disabled. Mitochondrial sequences were separated from the
nuclear assembly and maintained independently. Twenty-six short, high-depth
rDNA-associated contigs, totaling 82,735 bp, were archived separately as
unplaced repetitive sequences rather than classified as exogenous
contamination. The scaffold assembly was subsequently screened using NCBI
FCS-GX, and a 6,055-bp fish-derived scaffold was excluded after confirming
that it contained no annotated genes.

## Provenance

Procedure and rDNA treatment transcribed from:
`/Users/changyueyu/Documents/Codex/2026-08-24/ssolidus-single-cell-final-analysis/METHODS_current.md`,
section “Nuclear-genome assembly, polishing and curation”. Final release
statistics are also recorded in this repository's README.md.

The original Methods document remains as manuscript text; this file records
the assembly stage within the genome workflow. This is documentation, not
an executable reconstruction of the historical run.

## Runnable reconstruction

`01b_omnic_yahs_scaffolding.sh` supplies an executable reconstruction using
`yahs -q 10 -r 500000 --no-contig-ec`. Its BWA-MEM2 and duplicate-processing
choices are new implementation choices, not recovered historical parameters.
It supports paired Omni-C reads or an existing BAM and uses a new output
directory. It does not rerun contact-map curation or overwrite the final release.

Parameter references: https://github.com/c-zhou/yahs and
https://github.com/bwa-mem2/bwa-mem2 .
