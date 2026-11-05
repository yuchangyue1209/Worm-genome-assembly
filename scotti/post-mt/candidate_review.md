# BUSCO candidate investigation: evidence and caveats

This is a review set, not a rescue annotation and not a replacement for the solidus postprocessing workflow. Complete original inline scripts are not all duplicated here. Existing server TSV/FASTA/PAF outputs are authoritative for the recorded runs.

## odb10 comparison

Primary genome: 697 complete groups. All BRAKER proteins: 657 complete. Longest proteins: 656 complete.
98 genome-complete groups were not complete in all BRAKER proteins: 24 fragmented and 74 missing. Conversely, 58 complete groups were recovered by proteins from groups not complete in genome mode. Net difference: 40 complete groups.
The genome result reports 216/697 complete groups with internal stops. Simple exported-protein '*' counting reproduced only 29 and was rejected; no candidate classification was written from that failed audit.
The corrected BUSCO 6.0 audit used Miniprot ##PAF metadata, matched to exact selected protein IDs, and reproduced 697/216. Among the 98 candidates: 28 had stops in every selected copy; 70 had at least one selected copy without stops. No-stop status alone does not rule out frameshifts or incorrect models.
Output: `results/annotation_qc/03_busco_comparison/candidate_miniprot_stop_audit.tsv`. Script 21 reconstructs the validated procedure.

## Mapping to the current nuclear assembly

The 70 groups yielded 72 primary genomic regions, using BUSCO-reported locus bounds plus approximately 2 kb flanks. A candidate-region manifest preserves original coordinates. These windows are not curated CDS boundaries.
Initial Minimap2 asm5 output omitted -c. The reported 37 strong / 35 low-identity partition used approximate PAF match/length fields and was withdrawn. Do not use that classification.
Corrected mapping used -x asm5 -c --cs=short --secondary=yes -N 10. All 72 regions have strong mappings; no region has multiple strong mappings. The displayed lowest-identity hits were 100% coverage and 100% identity.
36 mapped windows overlap BRAKER mRNA spans; 36 do not. This is region-level overlap, not evidence of correct CDS reconstruction.
Files under `results/annotation_qc/04_busco_candidate_mapping`:
- primary_candidate_regions.fa and primary_candidate_regions.tsv
- candidates_vs_nuclear.baselevel.paf (use this corrected result)
- baselevel_braker_overlap.tsv
- unannotated_regions_evidence.tsv

## RNA and softmasking screen of 36 nonoverlapping regions

Unique primary alignment records were counted with samtools -F 2308 -q 255, using STAR's MAPQ 255 convention. Spliced records were identified by N in CIGAR. These are alignment records, not pairs or validated candidate introns.
35/36 windows had unique spliced records. candidate_042 had 6 unique records and 0 spliced records. Window softmasking ranged from approximately 19.5% to 83.2%.
Flanks or other overlapping transcription may supply these records; candidate strand, exact introns, complete CDS, and protein support have not been reconciled. Masking fraction alone does not establish why BRAKER omitted a model.
Decision: retain all evidence as an optional prioritized review set; do not automatically add models or restart training solely from these counts.

## Relationship to the solidus workflow

The local solidus script 07 used AGAT longest-isoform extraction, protein checks, and BUSCO metazoa_odb12. The README also records TSA/TagSeq refinement, addition of 187 supported genes, and removal of 62 unsupported extremely short models. Full construction code for all those refinements is not present in that directory.
For cotti, corresponding TSA/TagSeq evidence has not been established. The current main workflow proceeds through protein QC, Swiss-Prot, eggNOG, and InterProScan before conservative filtering. Do not describe this BUSCO investigation as the historical solidus gene-rescue method.
