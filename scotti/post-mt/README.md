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
|08_repeatmodeler.sh|Completed; do not launch a duplicate|
|09_raw_fastqc.sh|Raw RNA quality control; completed|
|10_sortmerna_pilot.sh|Reconstructed first-100,000-pair pilot; original pilot completed|
|11_run_rrna_chunks.py|Full paired screening with bounded chunks, validation, compression, checkpoints, and resume support|
|12_launch_full_rrna_filter.sh|Completed: 253 chunks; /work temporary storage|

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
The first failed launch left two empty output directories. Their emptiness was confirmed; the next supplied command removes only those empty directories with rmdir before launching the /work version. The full run subsequently completed successfully; see current_status.md.

The full script validates the original pilot, then filters sequential chunks without changing raw inputs. It deletes only its own successful chunk's temporary directory after validating and compressing the outputs. Failed chunks are retained.
It requires at least 100 GiB free before each chunk. Outputs go to `data/rnaseq/rrna_screened_v1/{rrna,non_rrna}`. Logs, counts, input records, and mate file lists go to `results/rnaseq_qc/05_sortmerna_full`.

Resume requires matching configuration: retain the two-million-pair chunk size rather than reverting to the ten-million default. Failed chunks with compressed outputs but no checkpoint require manual review; the script will not overwrite them.

## Subsequent scripts and current status

|Script|Purpose and status|
|---|---|
|13_repeatmasker.sh|Completed softmasking; 72.13% masked|
|14_filtered_fastqc_first_chunk.sh|First retained-chunk QC supplied; completion not reported|
|15_star_reference_and_index.sh|Completed nuclear-plus-mt index; omits sjdbOverhang at annotation-free indexing|
|16_star_rna_alignment.sh|Completed two-pass STAR mapping with sjdbOverhang 100 and XS tags|
|17_nuclear_bam.sh|Completed mt exclusion; 412 nuclear references and ptg000291l retained|
|18_braker3_etp_corrected.sh|Completed historical retry2 launch; do not rerun|
|18a_braker_initial_setup.sh|Historical initial configuration/protein copy; already done, do not repeat|

The earlier command included sjdbOverhang during annotation-free indexing and failed. Script 15 records the corrected successful method. RNA mapping uses the unmasked nuclear-plus-mt reference; BRAKER uses only the softmasked nuclear FASTA and nuclear BAM.
The BRAKER dependency check initially failed because its selected Python could not import Bio. The active BRAKER environment itself has Biopython 1.88. PATH shadowing by stringtie_env is the leading explanation; the corrected launch prioritizes braker3_annotation/bin and checks Bio after setting PATH. Retry2 subsequently completed successfully; see current_status.md.

See [current_status.md](current_status.md) for completed results, exact server paths, and the ptg000291l investigation. See [execution_record.md](execution_record.md) for earlier assembly decisions.

## Post-BRAKER scripts (updated 2026-11-04)

All scripts and comments are English. These reconstruct commands/results already discussed; updating local files does not run server jobs.

|Script|Purpose/status|
|---|---|
|19_protein_qc_and_longest.py|Completed protein QC and longest-isoform selection; requires a new --out directory|
|20_busco_proteins.sh|Completed raw10, longest10, and longest12 evaluations; existing BUSCO runs are protected|
|21_busco_comparison_and_stop_audit.py|Completed odb10 comparison and authoritative Miniprot stop audit; requires new --out|
|22_candidate_region_mapping.sh|Completed base-level mapping of previously extracted candidate regions; -c is essential|
|23_swissprot_diamond.sh|Completed all-isoform search using the solidus Swiss-Prot database|
|24_protein_support_audit.py|Completed short/long/X support audit; requires a new --out|
|25_eggnog_mapper.sh|Corrected installed CLI uses -i/-m; successful run not yet reported|
|26_interproscan_pending.sh|Guarded template only; no execution reported|

Example QC reconstruction in a new output directory:
```bash
python3 19_protein_qc_and_longest.py --out /work/cyu/annotation/scotti-annotation/results/annotation_qc/longest_reconstructed
```
Script 20 expects the original documented longest-protein path, not this reconstruction path. Existing completed results should be reviewed rather than rerun.
The candidate extraction, overlap, and regional RNA scripts originally supplied inline are summarized in candidate_review.md; they are not falsely represented as fully recovered server scripts here.

The primary-genome and cotti protein BUSCO comparison uses odb10 (954 groups). Comparison with the documented solidus final release uses odb12 (672 groups). Percentages from different lineages must not be treated as a change in annotation quality.
