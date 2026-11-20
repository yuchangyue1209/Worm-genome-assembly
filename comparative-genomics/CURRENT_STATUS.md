# Comparative genomics status — 2026-11-05

## Species design and existing work

The six-taxon core consists of S_solidus, S_cotti, L_intestinalis, D_latus, Spirometra_Aus1 and S_proliferum. Five Cyclophyllidea background taxa remain in the existing design: E_multilocularis, H_microstoma, E_granulosus, T_solium and T_multiceps. Do not interpret the core list as an instruction to delete these backgrounds.
The reported BUSCO batch contains ten species (all except cotti), each run with odb10 and odb12: 20 FINISH records. Actual scores and per-run inputs still need to be assembled from source summaries. Do not invent a score table from FINISH messages alone.
cotti has separate, validated candidate results and can be added explicitly for exploratory orthology, making eleven taxa. Final inference remains pending input curation and TE review.

## Completed S. proliferum data audit

SprJ_4.0, GCA_053814205.1 / PRJDB8966. Raw directory:
`/work/cyu/Schistocephalus_comparative_genomics/01_raw_data/S_proliferum_SprJ4/ncbi_dataset/data/GCA_053814205.1`.
29,231 genes; 29,231 mRNAs; 28,993 proteins; zero multi-mRNA genes. 238 mRNAs have no CDS; all other 28,993 have corresponding protein IDs. Zero FASTA proteins are unreferenced by CDS. These are audited no-CDS models, not automatically designated pseudogenes.
Audit file: `03_qc/S_proliferum_id_audit_v1/transcript_protein_audit.tsv`. Standardized raw links were created. Completion of the proposed primary-protein link should be checked on the server before launching new jobs.
The published TE-associated ID list and a reviewed non-TE set remain pending. The earlier chat cites 6,316 loci; their exact identities and mappings must be verified from supplementary source data, not inferred by subtracting that number or keyword filtering.

S_proliferum InterProScan was last observed running after approximately 37 minutes, with 24 CPUs, prefix `05_functional_annotation/interpro/S_proliferum`. Its actual log was not at the template path. No completion or match count has been reported. Find stdout via the live PID's /proc/<PID>/fd/1 if necessary; historical PIDs must not be reused blindly. Other S_proliferum functional/secretome completion states are unconfirmed in this chat.

## cotti handoff

Validated candidate: `/work/cyu/annotation/scotti-annotation/results/filtered_candidates/Scotti_filtered7_v1`.
9,234 genes / representatives; 11,464 isoforms and CDS sequences. GFF3 passed GenomeTools; odb12 C73.5% [S72.5%, D1.0%], F9.8%, M16.7%, n672, unchanged after seven exclusions.
Primary evidence: Swiss-Prot 6,545; eggNOG 7,351; member-database matches 8,490; integrated InterPro 7,499; any source 8,592. 642 representatives have no hits and are retained. Functional evidence is complete, but product naming, ncRNA work and final structural/TE review remain pending. Do not rerun searches merely to subset representative proteins.
Use scripts/00d_stage_scotti_candidate.sh to stage a separate candidate set. Stable 02_standardized_data remains unchanged. Candidate staged data are not a finalized release.

## BUSCO policy

Main cross-species reporting: the same local metazoa_odb12 release, same BUSCO version and protein mode, one protein per gene. Preserve odb10 results for historical/model investigations. Record dataset creation dates and software versions. Different dataset percentages are not directly comparable.
Script 02b now separates lineage-specific output and log directories, preventing an existing odb10 summary from causing an odb12 skip. Script 02d loops through both lineages. Script 02e discovers legacy and new summaries, preserves all source paths and flags duplicates rather than silently choosing a run. The reported batch is already completed: summarize before resubmitting.
The legacy 02c summary script remains available for a known single-lineage root; use 02e for the mixed historical tree. Its legacy quality labels are heuristics, not pass/fail biological validation.

## OrthoFinder boundary

Script 03 creates a new run directory, copies/freeze protein FASTAs, records source paths and hashes, and never reuses the old shared input folder. This avoids stale species symlinks. INCLUDE_SCOTTI_CANDIDATE=1 explicitly opts into a candidate run; results remain exploratory. Do not run CAFE5 directly from an unreviewed count matrix or interpret annotation absence as confirmed biological loss.
No new OrthoFinder, BUSCO or functional job was submitted by this organization task. Existing secretome/MEROPS scripts and parameters were preserved.

## Local verification of this organization update

All 54 shell/Python scripts across scotti/post-mt and comparative-genomics passed syntax checks, including embedded Python. Small synthetic fixtures verified mixed-lineage BUSCO collection, the 46-isoform/44-gene review split, and a simulated OrthoFinder input freeze that excludes stale files while preserving previous input directories. The functional-evidence integrator was previously tested with synthetic retained/excluded query IDs. These are code checks, not reruns of biological analyses. Server source files and active jobs were not accessed or modified.
