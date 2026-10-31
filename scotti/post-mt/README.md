# cotti workflow after mitochondrial assembly

Organized on 2026-10-28. Assembly root: `/work/cyu/assembly/worm/scotti`. Annotation root: `/work/cyu/annotation/scotti-annotation`.

These scripts were organized from commands and results supplied in the conversation. Steps without complete original commands are reconstructions, not copies of the original server scripts. No server analysis was rerun during this update.
Do not execute every numbered script in sequence: the directory includes completed experiments as well as the selected annotation workflow.

## Selected inputs

Nuclear annotation input: `final/Scotti_nuclear.mt_filtered.ragtag.provisional.fa`, containing 412 sequences and 991,317,046 bp.
This is the assembly after RagTag and before SAMBA, with four whole mitochondrial candidate sequences removed (89,757 bp). Large nuclear sequences with local mitochondrial matches were retained.

Validated HiFi reads: `reads/hifi/scotti_2.from_bam.fastq.gz`, containing 4,662,246 reads and 26,074,063,243 bp.
The older `scotti_2.fastq` is malformed and must not be used by these scripts. Whether its corruption arose during merging or transfer remains unresolved.

Mitochondrial working package: `mitochondria/08_working_release/Scotti_mt_working_v1`. Assembly length: 15,437 bp. Annotation: 12 protein-coding genes, 22 tRNAs, and 2 rRNAs. The GFF3 passed GenomeTools validation but remains a working candidate.

## Script index

|Script|Purpose and status|
|---|---|
|00_config.sh|Shared server paths; does not launch analysis|
|01_samba_hifi_experiment.sh|Reconstructed HiFi experiment using `-d asm -m 5000`; experiment completed, output not selected for annotation|
|02_compare_assemblies.py|Assembly statistics and exact sequence comparisons independent of names, order, and orientation; exports changed sequences|
|03_primary_depth_assessment.sh|Reconstructed primary-contig PAF/pbcstat assessment; assessment completed|
|04_primary_busco.sh|Reconstructed primary BUSCO command; assessment completed|
|05_map_hifi_to_samba.sh|HiFi mapping to new SAMBA output; historical experiment completed|
|06_mt_source_tracking.sh|BLAST evidence across primary, RagTag, and SAMBA assemblies|
|07_separate_reviewed_mt.py|Reconstructs reviewed RagTag-only mitochondrial separation in a new output directory|
|08_repeatmodeler.sh|Started on the server; completion not yet reported; do not launch a duplicate|
|09_raw_fastqc.sh|Raw RNA quality control; completed|
|10_sortmerna_pilot.sh|Reconstructed first-100,000-pair pilot; original pilot completed|
|11_run_rrna_chunks.py|Full paired screening with bounded chunks, validation, compression, checkpoints, and resume support|
|12_launch_full_rrna_filter.sh|Uses two million pairs per chunk and /work temporary storage; successful full-run start not yet confirmed|

Run shell scripts with `bash`, not `sh`. Activate the environment stated in each script first.
Directories ending in `_reconstructed` protect historical results. Scripts refuse existing output directories by default.
The comparison and separation Python scripts require Biopython. The full RNA filtering script requires Python's standard library and an installed SortMeRNA executable.

Comparison example, on the server in the genome_assembly environment:

```bash
python3 02_compare_assemblies.py /work/cyu/assembly/worm/scotti/qc/comparison_reconstructed \
 before=/work/cyu/assembly/worm/scotti/raw/scotti_2_ragtag_500min_direct.fasta \
 new=/work/cyu/assembly/worm/scotti/scaffolding/samba_hifi_asm_m5000/input.fa.scaffolds.fa
```

## RNA screening and storage failure

Kyle identified the raw input as uncleaned paired RNA-seq: 505,478,106 reads per mate, 101 bp each.
Quality and adapter checks passed. GC, duplication, and overrepresented-sequence results warranted investigation; most overrepresented sequences matched solidus nuclear rRNA. High duplication alone is not a reason to deduplicate RNA-seq.

Reference: `results/rnaseq_qc/02_overrepresented_check/solidus_nuclear_rrna.fa`, containing 761 solidus annotation intervals. This cross-species reference cannot guarantee complete removal of cotti rRNA or distinguish host-derived rRNA.

Pilot validation confirmed mate identities and conservation of all 100,000 input pairs:

- rRNA-associated output: 74,077 pairs (74.08%).
- Retained output: 25,923 pairs (25.92%).

The 73.48% in the SortMeRNA log is a per-read hit rate. With `--paired_in`, a hit in either mate sends both mates to the rRNA-associated output. Retained reads should not be described as pure mRNA.

Reading `/scratch/cyu` and writing a temporary file both failed with EIO. Free space reported by df does not establish that the filesystem is usable. An administrator should investigate; these scripts do not attempt filesystem repairs.
The first failed launch left two empty output directories. Their emptiness was confirmed; the next supplied command removes only those empty directories with rmdir before launching the /work version. Full-run success has not yet been reported.

The full script validates the original pilot, then filters sequential chunks without changing raw inputs. It deletes only its own successful chunk's temporary directory after validating and compressing the outputs. Failed chunks are retained.
It requires at least 100 GiB free before each chunk. Outputs go to `data/rnaseq/rrna_screened_v1/{rrna,non_rrna}`. Logs, counts, input records, and mate file lists go to `results/rnaseq_qc/05_sortmerna_full`.

Resume requires matching configuration: retain the two-million-pair chunk size rather than reverting to the ten-million default. Failed chunks with compressed outputs but no checkpoint require manual review; the script will not overwrite them.

## Pending dependencies

RepeatModeler completion → inspect the repeat library → RepeatMasker softmasking.
Full RNA screening → post-filter QC → STAR splice-aware mapping → RNA/protein evidence for BRAKER.
For 101 bp reads, STAR sjdbOverhang should be 100; do not copy the solidus value of 149.
STAR, RepeatMasker, and BRAKER are not recorded as executed. Their inputs and environments must be verified after the preceding steps complete.

See [execution_record.md](execution_record.md) for results and decisions.
