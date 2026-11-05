# S. cotti mitochondrial assembly and annotation

Compiled on 2026-10-27 for **S. cotti / scotti_2 / PacBio HiFi**.
The original `workflow-scripts/02_mitochondrial_assembly_annotation.sh` describes
**S. solidus / ONT** and remains unchanged.
This directory contains reproducible workflow code and results reported in the
analysis session. Large data files and actual analysis outputs remain on stickleback.
Local validation covered syntax and small logic checks; assembly and annotation
were not rerun locally.

## Current status

- The BAM-derived FASTQ passed a complete format check. The old `scotti_2.fastq` is malformed and must not be reused.
- A single 15,437-bp circular mt candidate was assembled, with read support for the junction and repeat region.
- MITOS2 and ARAGORN have been run. Homology supports adding rrnL and trnH and revising the cox3 boundaries.
- **A revised GFF has not been generated, validated, or released.** `curation_candidates.tsv` records proposed changes, not final annotations.
- OH predictions are not confirmed D-loops. Do not add genomic A residues to the FASTA to complete an inferred stop codon.
- These annotation issues do not, by themselves, justify rerunning nuclear hifiasm assembly.

## Data paths and sample identity

```text
/work/cyu/assembly/worm/scotti/
  raw/scotti_2_ragtag_500min_direct.fasta.scaffolds.simpleheaders.fa
  reads/hifi/scotti_2_merged.bam              # Two sequencing runs of the same scotti_2 individual
  reads/hifi/scotti_2.from_bam.fastq.gz       # Validated replacement FASTQ
  reads/hifi/scotti_2.fastq                  # Malformed original; retained for provenance only
  reads/hifi/m84122_240202_182632_s2.hifi_reads.bam  # scotti_1; do not combine with scotti_2
  excluded/scotti_mt_candidates.fa
  final/Scotti_nuclear.mt_filtered.provisional.fa
  mitochondria/
    01_read_extraction_from_bam/
    03_flye_hifi_random700/assembly/assembly.fasta
    04_validation/
    05_mitos2/result.gff
```

The project's own solidus assembly is
`/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa` (15,436 bp).
The published reference is MW602521.1 (13,947 bp), originally stored as
`/work/cyu/assembly/worm/mitochondria/Ssolidus_mt_ref.fa`.
These are distinct references. The 1-bp total-length difference between the
project's solidus and cotti assemblies does not establish sequence identity or
equal repeat copy number.

## Running the workflow

Run these scripts in **Bash on stickleback**, not on the local Mac. Copy this
whole directory to the server, then select the required stage inside tmux.
Do not automatically repeat completed stages.
There is no master script that reruns the entire workflow. Major computational
stages require a fresh output directory and stop if it already exists.
Comparison stages regenerate the documented derived alignment files without
modifying the input FASTA or GFF.

```bash
conda activate genome_assembly
# Change to this directory on the server, then select the required stage:
bash scripts/01_bam_to_fastq.sh
# bash scripts/02_mt_candidates.sh
# bash scripts/03_extract_mt_reads.sh
# bash scripts/04_flye_random700.sh
# bash scripts/05_validate_assembly.sh

# Annotation and gene-boundary checks use a separate environment:
# conda activate mitos2
# bash scripts/06_mitos_aragorn.sh
# bash scripts/07_compare_annotations.sh
```

| Script | Purpose and notes |
|---|---|
| `00_config.sh` | Shared paths; sourced by other scripts; ROOT can be overridden |
| `01_bam_to_fastq.sh` | Unmapped BAM quickcheck -u, complete export, and full FASTQ validation |
| `02_mt_candidates.sh` | Initial BLAST, interval-union coverage, provisional separation of four reviewed candidates, and self-alignment; reconstruction outputs go to qc/mt_reproduced |
| `03_extract_mt_reads.sh` | map-hifi; retain mapped primary reads without a MAPQ threshold |
| `04_flye_random700.sh` | Randomly select 700 reads with seed 42; successful Flye HiFi route |
| `05_validate_assembly.sh` | Align both solidus references; rotate the sequence; assess junction and repeat-region read support |
| `06_mitos_aragorn.sh` | runmitos + refseq89m + table 9 + --best; ARAGORN -mt -gc9 |
| `07_compare_annotations.sh` | Extract rrnL/trnH/cox3 from the project's solidus annotation; retrieve the published GenBank record; translate alternative cox3 boundaries |
| `mt_helpers.py` | Supporting data processing; annotation functions require Biopython |

MITOS2 configuration used in this project:

```bash
conda activate mitos2
# If the conda-protect plugin blocks an operation, temporarily use:
# CONDA_NO_PLUGINS=true conda activate /home/cyu/.conda/envs/mitos2
command -v runmitos
# /home/cyu/.conda/envs/mitos2/bin/runmitos
# Database base directory: /work/cyu/assembly/worm/mitochondria/mitos_db
# Database version directory: refseq89m
```

Use `runmitos`, not the unavailable `runmitos.py`. Do not reuse the solidus
ONT/Medaka parameters for cotti HiFi reads.
Commands and checks have been organized into reusable scripts rather than
preserved as a verbatim terminal transcript. Failed attempts are documented in
`history_and_caveats.md`.

## Reported results

These values were reported from server runs, not recalculated locally.

### Reads and candidate sequences

- BAM export: 4,662,246 reads; 26,074,063,243 bp; mean length 5,592.6 bp; longest read 75,178 bp.
- Four candidates: scaffold206/312/349/410, totaling 89,757 bp; none contains Ns.
- Original assembly: 413 sequences, 991,407,400 bp. Provisional nuclear assembly: 409 sequences, 991,317,643 bp.
- These IDs apply only to this particular simpleheaders assembly. Mt-like segments within long nuclear scaffolds were not directly deleted.
- Extracted candidate set: 11,127 reads, 49,450,194 bp; NUMT-derived reads may remain.

### Flye assembly and read support

- Assembly of all candidate reads failed: overlap-based coverage 2121; zero disjointigs.
- Random 700-read subset (seed 42) succeeded: contig_1, 15,437 bp, coverage 211, circ. Y, repeat N, mult. 1.
- After rotation, the original junction is at the 0-based boundary 7719. **1,851 reads** span at least 1 kb on each side with MAPQ >=20 and alignment-wide identity >=95%.
- The candidate repeat region is approximately 8.3–10.3 kb in the original coordinates. **602 qualifying reads** span 7300–11300.
- Of these, 599 have no insertion/deletion >=100 bp in the tested region; three have deletions (525 bp twice, 526 bp once).
- The 599/602 proportion is conditional on the alignment filters, not an estimate of sample-wide heteroplasmy. Divergent and split alignments may be excluded.
- Junction support alone does not independently establish sequence origin or every repeat copy. Annotation and final review remain necessary.

### Annotation and boundary review

MITOS2 reported `missing: atp8 trnH rrnL OL` and `duplicated: 2x OH`.
ARAGORN reported 27 tRNA candidates, not 27 confirmed genes.
The ARAGORN trnH candidate at 5364–5428 on the minus strand overlaps the original
MITOS2 cox3 prediction at 5366–6016 by 63 bp.
Broad OH predictions and repeated hits are not confirmed replication origins
and must not be described directly as two D-loops.

Published MW602521.1 annotations, converted to 1-based inclusive coordinates:

- cox3: 7200–7842, plus strand, 643 bp; terminal T, annotated as completing TAA through addition of A residues.
- trnH: 7843–7906, plus strand, 64 bp; adjacent to cox3.

Supported cotti revision candidates:

| Gene | cotti coordinates, 1-based inclusive | Strand | Evidence |
|---|---|---|---|
| rrnL | 12124–13172 | − | Full query coverage of the project's solidus rrnL; 90.190% identity |
| trnH | 5365–5427 | − | Full published-reference query coverage, 96.875% identity; ARAGORN boundaries differ by approximately 1 bp at each end |
| cox3 | 5428–6070 | − | Full published-reference query coverage, 81.493% identity; GTG start, 214 complete codons, terminal T, no internal stop |

Direct translation of the proposed cox3 sequence begins with V; if confirmed as
the initiation codon, the protein begins with M.
The original MITOS2 cox3 prediction has ATG/TAA, but a complete ORF alone does not
establish correct boundaries. Published-reference homology supports the revised range.
The proposed changes have not been written into a GFF or subjected to submission
validation. The rrnL homology-based boundaries, tRNA boundaries, and complete
gene set still require final review.
A revised GFF must preserve consistent gene/transcript/exon/CDS relationships,
Parent/ID attributes, phases, and incomplete-stop annotations; changing one
coordinate row alone is insufficient.

## Remaining work

1. Record the reviewed annotation changes in a new file with provenance; preserve the original MITOS2 GFF.
2. Review the complete gene set, reading frames, rRNAs, tRNAs, and overlaps. Do not promote every ARAGORN prediction to the final annotation.
3. Validate feature hierarchy and coordinates before generating the final mt release.
4. Handle SAMBA settings and nuclear assembly validation in the separate nuclear workflow. Mt annotation uncertainty alone does not imply that nuclear reassembly is needed.

## References

- Flye: https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md
- minimap2: https://github.com/lh3/minimap2
- samtools quickcheck: https://www.htslib.org/doc/samtools-quickcheck.html
- Genetic code 9: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
- Published mitochondrial record: https://www.ncbi.nlm.nih.gov/nuccore/MW602521.1
