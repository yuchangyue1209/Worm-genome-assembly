# Current annotation status

Updated from user-supplied server output through 2026-11-04 post-BRAKER QC and the eggNOG CLI correction. Server logs label these runs Oct 31, 2026; this record preserves that label without changing it to the local client date. Local organization does not rerun server analyses.

## Completed repeats

RepeatModeler completed in 18:09:14. Preserve Scotti_ragtag_nuclear-families.fa, Scotti_ragtag_nuclear-families.stk, and Scotti_ragtag_nuclear-rmod.log under `repeats/01_repeatmodeler_fixed`.
RepeatMasker 4.2.4 with rmblastn 2.17.1+ completed against the cotti repeat library.
Input: 412 sequences, 991,317,046 bp, GC 44.26%. Masked: 715,044,557 bp (72.13%). LINEs: 45.33%; unclassified: 17.51%; total interspersed repeats: 69.93%.
BRAKER genome: `/work/cyu/assembly/worm/scotti/repeats/02_repeatmasker/Scotti_nuclear.mt_filtered.ragtag.provisional.fa.masked`.

## Completed full RNA screening

253 chunks, 505,478,106 input pairs. rRNA-associated: 378,734,366 pairs (74.93%). Retained: 126,743,740 pairs (25.07%). All chunks reported PASS and the run reported COMPLETE.
Outputs: `/work/cyu/annotation/scotti-annotation/data/rnaseq/rrna_screened_v1`.
Reports and mate input lists: `/work/cyu/annotation/scotti-annotation/results/rnaseq_qc/05_sortmerna_full`.
Pair policy: if either mate hits the solidus-derived reference, the pair goes to the rRNA-associated group. These categories do not establish a pure mRNA library or the original library preparation protocol.
Post-filter FastQC was supplied for the first chunk; its results have not been provided. Trimming has not been reported.

## Completed STAR mapping

STAR 2.7.11b. Reference combines the unmasked RagTag-only nuclear FASTA with the independently assembled mitochondrial FASTA, renamed Scotti_mt: 413 sequences, 991,332,483 bp.
Reference: `results/star/01_index/reference.nuclear_plus_mt.fa`.
Successful index: `results/star/01_index_no_gtf`. Initial annotation-free indexing with sjdbOverhang 100 failed; corrected indexing omitted that parameter.
Mapping: two-pass Basic, sjdbOverhang 100, outSAMstrandField intronMotif, coordinate-sorted BAM.
Input: 126,743,740 pairs; average combined read length 202 bp.
Unique: 77,116,760 pairs (60.84%); multiple loci: 12,646,541 (9.98%); total accepted mapping: 70.82%.
Too many loci: 1.69%; unmapped too short: 23.46%; unmapped other: 4.02%. “Too short” is STAR alignment scoring/length terminology, not proof of short input reads.
Mismatch rate: 0.29%. Average mapped length: 201.27 bp.
Splice counts: total 17,562,594; sjdb 17,393,434; GT/AG 17,222,783; GC/AG 315,626; AT/AC 5,727; noncanonical 18,458.
The sjdb category includes first-pass junctions and is not independent annotation validation.
Combined BAM: `/work/cyu/annotation/scotti-annotation/results/star/02_rna_alignment/Aligned.sortedByCoord.out.bam`.
Scotti_mt idxstats: 26,200,024 mapped alignment records; this is not a pair count.

## Completed nuclear BAM and XS checks

BRAKER BAM: `/work/cyu/annotation/scotti-annotation/results/star/02_rna_alignment/Scotti_nuclear.rnaseq.bam`.
Removed mapped mt records, unmapped records, and nuclear records whose mate reference is mt. Removed only Scotti_mt from the sequence dictionary. Retained ptg000291l and all 412 nuclear sequence names. Retained nuclear secondary alignments.
Records: total 217,684,946; primary 156,963,786; secondary 60,721,160; supplementary 0. Primary records represent 78,481,893 pairs.
All retained primary records were reported mapped and properly paired. These are post-filter statistics, not original-library mapping rates.
Spliced records: 26,379,328; records with XS: 26,379,328. XS check passed.

## ptg000291l investigation: retained, pending annotation review

Contig length 3,243,445 bp. idxstats initially showed 13,847,428 mapped records, including secondary records.
Primary mapped records: 4,225,130. NH=1: 2,436,608; NH>1: 1,788,522.
BLAST to solidus-derived rRNA intervals, identity >=80% and alignment length >=100 bp: 48 merged intervals covering 19,553 bp (0.60%). Only 1,753 primary records overlapped these intervals (0.04%). This does not support explaining the whole contig's high signal as identified rRNA; unrecognized rRNA remains possible.
10 kb alignment-start bins: 2,610,000–2,620,000 had 2,594,355 records; 2,630,000–2,640,000 had 794,763. Together about 80.2% of primary starts. These are start counts, not per-base depth, and whole-contig NH fractions do not establish hotspot-specific fractions.
RepeatMasker reported 84 overlapping annotations in the 2,600,001–2,670,000 window, including LINE/CR1, LINE/RTE-BovB, LTR/Gypsy, and other families.
Repeated family patterns suggest a candidate approximately 13–14 kb structural repeat. SINE/Alu labels do not establish human contamination.
Self-alignment asm20 -DP: one independent reciprocal match of about 8.5 kb, 91.3% identity, offset 26,857 bp. This supports local duplication, not full validation of four tandem units or evidence of assembly error.
Decision: retain sequence and RNA evidence, use softmasking, and review predictions in 2.60–2.67 Mb after annotation. No whole-contig exclusion, reassembly, or full RNA refiltering was performed on this basis.

## BRAKER setup and completed run

Core method matches solidus: BRAKER3 ETP, softmasked nuclear FASTA, RNA BAM with XS, and the same filtered 128,216-sequence protein evidence library.
Protein source: `/work/cyu/annotation/ssol-annotation/results/evidence/braker_proteins.filtered.fa`.
Copied input: `/work/cyu/annotation/scotti-annotation/results/braker3_etp_v1/input.proteins.fa`.
Private AUGUSTUS configuration: `/work/cyu/annotation/scotti-annotation/config/augustus_scotti_v1`.
Training species: scotti_cyu_etp_v1; threads: 32.
GeneMark: `/work/cyu/annotation/software/GeneMark-ETP/bin`.
ProtHint: `/work/cyu/annotation/software/GeneMark-ETP/bin/gmes/ProtHint/bin`.
StringTie 2.1.7: `/home/cyu/.conda/envs/stringtie_env/bin/stringtie`.
BRAKER environment Python: `/home/cyu/.conda/envs/braker3_annotation/bin/python3`, Bio version 1.88 verified.

Initial run braker3_etp_v1 and retry1 stopped at the Bio import dependency check. ProtHint SyntaxWarnings were not the reported fatal error. The earlier PATH prepended stringtie_env/bin, potentially selecting a different Python; this is the leading diagnosis, not directly confirmed by a captured failing interpreter path.
Corrected script 18 places braker3_annotation/bin first, then stringtie_env/bin, and validates Python after setting final PATH.
Successful run: `/work/cyu/annotation/scotti-annotation/results/braker3_etp_v1_retry2`.
**Retry2 completed successfully.** The log ends with BRAKER RUN FINISHED. Final output files were reported at Oct 31 10:15 after an Oct 31 02:12:56 start. Preserve failed logs; do not rerun the completed job. Raw predictions: 9,241 genes, 11,471 mRNAs, 85,079 exons and 85,079 CDS features. Both protein and coding-sequence FASTAs contain 11,471 sequences. These are pre-filter predictions, not a final release.

## Post-BRAKER QC (reported 2026-11-04)

|Protein set / lineage|Complete|Single|Duplicated|Fragmented|Missing|n|
|---|---:|---:|---:|---:|---:|---:|
|All isoforms / metazoa_odb10|68.9%|58.6%|10.3%|5.7%|25.5%|954|
|Longest / metazoa_odb10|68.8%|66.8%|2.0%|5.7%|25.6%|954|
|Longest / metazoa_odb12|73.5%|72.5%|1.0%|9.8%|16.7%|672|

Longest-isoform selection used exact GFF3 mRNA ID/Parent mapping to braker.aa and retained one longest protein per gene (ID tie-break). It is not guaranteed identical to AGAT's selection in solidus. Original isoforms remain intact.
Longest set: 9,241 sequences, 0 internal stops, 9 containing X, 44 shorter than 50 aa, 452 shorter than 100 aa, median 403 aa, maximum 8,418 aa.
Solidus final revised longest proteome: odb12 C80.4%, F7.0%, M12.6%, n672. The two species' results are at different annotation stages; no direct causal conclusion follows from the difference.

## Swiss-Prot and protein support

All 11,471 BRAKER proteins were searched with DIAMOND, sensitive, E-value 1e-5, max-target-seqs 5, using the existing solidus Swiss-Prot database. Completed in about 8.9 seconds; 38,695 reported alignments; 8,206 proteins with hits (71.54%).
Results: `results/functional/01_swissprot/scotti_vs_swissprot.{top5,best_hit}.tsv`.
Audit: `results/annotation_qc/06_protein_support_audit/protein_support.tsv`.

|Category, all isoforms|Total|Swiss-Prot hit|
|---|---:|---:|
|<50 aa|46|2|
|>5,000 aa|34|28|
|Contains X|9|6|

Hits do not establish full-length support, and absence of a Swiss-Prot hit does not establish a false gene. No proteins have been removed. Short/long/X models require additional evidence.

## Functional annotation pending

The first eggNOG launch failed at argument parsing: this installed emapper.py rejects --input and --method. Corrected command uses -i and -m diamond, with --data_dir /work/cyu/databases/eggnog_mapper, all 11,471 proteins, prefix scotti_braker, output results/functional/02_eggnog, and 32 CPUs. A corrected run command was supplied; successful start/completion has not yet been reported.
Do not use --override or copy solidus's expected row count of 8,177 into cotti validation. Script 25 preserves existing logs and refuses existing emapper result files. Inspect any active job before launching.
InterProScan remains pending. Script 26 is a guarded template adapted from the solidus installation; verify its Java/runtime requirements before execution. It removes only a terminal stop from the search copy, preserves IDs and X, and rejects unsupported residues. Original proteins remain unchanged.

## Current next steps

Confirm corrected eggNOG execution; then validate/launch InterProScan when resources permit. Integrate evidence and review unsupported short models and potential long-gene fusions before any conservative filtering. The 70 BUSCO candidate groups remain a separate review set, not automatically rescued gene models. See candidate_review.md.
