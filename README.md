# *Schistocephalus solidus* Genome Assembly and Annotation

This repository records the workflow used to generate the 2026 nuclear and
mitochondrial genome assemblies and the evidence-refined genome annotation of
*Schistocephalus solidus*. It documents the completed analysis, release-quality
control, and the remaining steps required for a public functional annotation.

## Workflow overview

1. Filter Oxford Nanopore reads and assemble the nuclear genome with hifiasm.
2. Polish with Racon and Medaka and remove redundant haplotigs with purge_dups.
3. Assemble and annotate the mitochondrial genome separately.
4. Build a species-specific repeat library and soft-mask the nuclear assembly.
5. Align public RNA-seq data and generate intron hints.
6. Run BRAKER3 in ETP mode with RNA and protein evidence.
7. refine gene models with TSA and TagSeq evidence.
8. perform contamination, protein, rDNA, ncRNA, and BUSCO quality control.
9. annotate proteins with Swiss-Prot, eggNOG-mapper, and InterProScan.
10. integrate conservative product names and prepare the NCBI release.

## Numbered workflow scripts

| Script | Purpose |
|---|---|
| `01_ont_nuclear_assembly.sh` | ONT filtering, hifiasm assembly, Racon/Medaka polishing, purge_dups, and assembly QC records |
| `02_mitochondrial_assembly_annotation.sh` | Mitochondrial read extraction, Flye assembly, Medaka polishing, MITOS2, and ARAGORN |
| `03_document_mt_removal.sh` | Document the separation of mitochondrial and nuclear scaffolds |
| `04_repeat_annotation.sh` | RepeatModeler2 discovery and RepeatMasker soft masking |
| `05_rnaseq_evidence.sh` | STAR alignment of 15 PRJNA304161 libraries and generation of intron hints |
| `06_braker3_etp.sh` | BRAKER3 ETP structural annotation using RNA and protein evidence |
| `07_annotation_postprocess_qc.sh` | Canonical-isoform extraction, protein validation, and BUSCO |
| `08_structural_release_manifest.sh` | Validate the evidence-refined structural annotation release |
| `09_fcs_gx_contamination_audit.sh` | Document FCS-GX screening and removal of the fish contaminant scaffold |
| `10_final_protein_and_swissprot_audit.sh` | Audit unusual proteins and summarize Swiss-Prot coverage |
| `11_swissprot_diamond.sh` | Run DIAMOND against UniProtKB/Swiss-Prot and select best hits |
| `12_ncrna_primary_scans.sh` | Record tRNAscan-SE, barrnap, and whole-genome Rfam scans |
| `13_ncrna_annotation_and_audit.sh` | Build the conservative ncRNA release and audit CDS overlaps |
| `14_release_v1_manifest.sh` | Validate current assembly and annotation components and write checksums |

The first ONT script preserves exploratory command branches and should be read
as a historical workflow record. Expensive primary scans should not be repeated
unless their input sequences or database versions change.

## Current nuclear assembly

Final FASTA:

```text
/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
```

| Metric | Value |
|---|---:|
| Assembly size | 846,268,314 bp |
| Scaffolds | 166 |
| Scaffold N50 | 99,718,931 bp |
| Longest scaffold | 154,525,006 bp |
| GC content | 44.07% |
| Ambiguous bases | 5,122 bp |
| Repetitive sequence | 68.64% |

Current SHA-256:

```text
16455061212d97b94e3f20abcedee353719dc651026adb39744cdb577ca3b2f6
```

FCS-GX identified one 6,055-bp scaffold (`scaffold_61`) as fish-derived,
assigned to *Esox lucius*, and recommended `EXCLUDE`. The scaffold contained no
records in the final GFF3 or GTF. Its removal reduced the assembly from 167 to
166 scaffolds without changing the protein-coding annotation or proteome.

## Mitochondrial assembly

```text
/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa
/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.gff3
```

The mitochondrial genome was assembled and annotated separately from the
nuclear genome. Its annotation includes the 12 protein-coding genes expected
for cestode mitochondrial genomes.

## Repeat annotation

RepeatModeler2 and RepeatMasker estimated that 68.64% of the nuclear assembly is
repetitive. This value describes genome composition; it should not be treated
as evidence that repeat transcription alone explains read-mapping behavior.

## Structural annotation

BRAKER3 was run in ETP mode with public RNA-seq, TSA proteins, and proteins from
related cestodes. The recovered BRAKER3 intermediate contained 9,671 genes and
12,724 transcripts. Its canonical proteome BUSCO result was:

```text
C:79.6%[S:79.2%,D:0.4%],F:6.1%,M:14.3%,n:672
```

TSA evidence corrected 3,707 UTRs, 187 high-confidence new genes were added,
and TagSeq evidence corrected 4,340 3-prime ends.

Final nuclear annotation files:

```text
/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3
/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.checked.gtf
/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_final.proteins.fa
```

| Feature | Count |
|---|---:|
| Protein-coding gene loci | 9,858 |
| Transcripts/protein isoforms | 12,911 |
| Exons | 105,217 |
| CDS features | 104,725 |
| 5-prime UTR features | 4,247 |
| 3-prime UTR features | 7,830 |

Final proteome BUSCO (`metazoa_odb12`):

```text
C:80.4%[S:79.9%,D:0.4%],F:7.0%,M:12.6%
```

## Protein audit

The final proteome contains 12,911 sequences totaling 7,894,936 amino acids,
with a mean length of 611.5 aa and a maximum length of 8,390 aa. No internal
stop codons were detected. The audit identified 77 proteins shorter than 50 aa
and 29 proteins longer than 5,000 aa. One short and 26 long proteins have
Swiss-Prot hits. Among the remaining short models, 14 have TagSeq support and
62 remain low-confidence review candidates. No model should be removed solely
because of length.

## rDNA and nuclear ncRNA annotation

Twenty-six short rDNA-associated contigs totaling 82,735 bp were excluded during
assembly cleanup. The 28,068-bp `ptg001004l` contig had 18.31-fold median depth
and contained redundant LSU/5.8S/SSU rDNA units. It was not restored because
the retained assembly contains complete high-copy rDNA arrays, including about
19 units in the first 250 kb of `scaffold_9`.

The conservative nuclear ncRNA release is:

```text
/work/cyu/annotation/ssol-annotation/results/ncrna_final/release/Ssolidus_nuclear.ncRNA.v1.0.gff3
```

| ncRNA type | Count |
|---|---:|
| Near-complete rRNA candidates | 761 |
| Conservative tRNAs | 352 |
| snRNAs | 32 |
| snoRNAs | 2 |
| miRNA | 1 |
| RNase P RNA | 1 |
| Total | 1,149 |

The high rRNA count reflects tandem high-copy arrays and does not represent 761
independent functional genes. A broader 2,612-feature Rfam candidate set is
retained as supplementary annotation. No conservative ncRNA feature overlaps a
protein-coding CDS on either strand.

## Functional annotation

### Swiss-Prot

DIAMOND against UniProtKB/Swiss-Prot release 2026_02 produced significant hits
for 9,172 of 12,911 proteins (71.04%).

### eggNOG-mapper

eggNOG-mapper 2.1.13 with eggNOG database 5.0.2 annotated 8,177 proteins:

| Annotation field | Proteins |
|---|---:|
| Any eggNOG annotation | 8,177 |
| Description | 8,103 |
| Preferred name | 6,877 |
| GO terms | 6,896 |
| EC number | 2,267 |
| KEGG KO | 6,554 |
| Pfam | 7,947 |

Swiss-Prot and eggNOG jointly cover 9,538 unique proteins (73.88%).

### InterProScan

InterProScan 5.78-109.0 was installed and validated with a 10-protein test. The
full 12,911-protein analysis is the current functional-annotation step. Its
results will be integrated with Swiss-Prot and eggNOG after completion.

## Single-cell validation note

The nuclear annotation and separate mitochondrial annotation were combined to
build a 10x reference as an independent validation of annotation usability.
Including mitochondrial sequence substantially improved transcriptome mapping,
showing that the previously excluded mitochondrial genome was an important
source of unmapped reads. These single-cell results are validation evidence and
are not the primary focus of this assembly repository.

## Remaining release steps

1. Complete the full InterProScan analysis.
2. Integrate Swiss-Prot, eggNOG, and InterPro evidence.
3. Review the 62 short low-confidence protein models.
4. Assign conservative product names.
5. Merge protein-coding and conservative ncRNA annotations.
6. assign stable `locus_tag`, `transcript_id`, and `protein_id` values.
7. run table2asn and NCBI validation.
8. freeze release v1.0 and generate checksums.
