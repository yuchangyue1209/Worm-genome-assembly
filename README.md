# *Schistocephalus solidus* Genome Assembly and Annotation

This repository records the reproducible workflow used to produce the 2026
*Schistocephalus solidus* nuclear and mitochondrial assemblies, the
evidence-refined protein-coding annotation, the conservative ncRNA annotation,
and the integrated functional annotation.

The scientific release is complete. Registration of an NCBI locus-tag prefix,
submission-file generation, and NCBI validation are intentionally treated as a
separate submission stage.

## Workflow scripts

| Script | Purpose |
|---|---|
| `01_ont_nuclear_assembly.sh` | ONT read filtering, hifiasm assembly, Racon/Medaka polishing, purge_dups, and assembly QC |
| `02_mitochondrial_assembly_annotation.sh` | Mitochondrial read extraction, Flye/Medaka assembly, MITOS2, and ARAGORN |
| `03_document_mt_removal.sh` | Document separation of mitochondrial and nuclear sequence |
| `04_repeat_annotation.sh` | RepeatModeler2 discovery and RepeatMasker soft masking |
| `05_rnaseq_evidence.sh` | RNA-seq alignment and intron-hint generation |
| `06_braker3_etp.sh` | BRAKER3 ETP structural annotation |
| `07_annotation_postprocess_qc.sh` | Isoform extraction, protein validation, and BUSCO |
| `08_structural_release_manifest.sh` | Validate the pre-functional structural release |
| `09_fcs_gx_contamination_audit.sh` | FCS-GX screening and fish-contaminant exclusion audit |
| `10_final_protein_and_swissprot_audit.sh` | Short/long protein audit and Swiss-Prot summary |
| `11_swissprot_diamond.sh` | DIAMOND search against UniProtKB/Swiss-Prot |
| `12_ncrna_primary_scans.sh` | Optional tRNAscan-SE, barrnap, and whole-genome Rfam scans |
| `13_ncrna_annotation_and_audit.sh` | Conservative ncRNA construction and CDS-overlap audit |
| `14_release_v1_manifest.sh` | Validate the completed filtered structural and functional release |
| `15_eggnog_mapper.sh` | Run or validate eggNOG-mapper annotation |
| `16_interproscan.sh` | Run or validate the full InterProScan annotation |
| `17_integrate_functional_annotations.sh` | Integrate Swiss-Prot, eggNOG, and InterPro evidence and product names |
| `18_merge_functional_ncrna_release.sh` | Merge protein-coding and ncRNA GFF3 and sanitize structural fields |
| `19_validate_final_release.sh` | Perform final coordinate, sequence-ID, feature-count, and checksum validation |
| `20_publish_shared_release.sh` | Build the versioned shared pre-NCBI scientific release package |

Expensive primary searches are opt-in and should not be repeated unless their
input sequence or database version changes.

## Final nuclear assembly

```text
/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
```

| Metric | Final value |
|---|---:|
| Assembly size | 846,268,314 bp |
| Nuclear scaffolds | 166 |
| Scaffold N50 | 99,718,931 bp |
| Longest scaffold | 154,525,006 bp |
| GC content | 44.07% |
| Ambiguous bases | 5,122 bp |
| Repetitive sequence | 68.64% |

FCS-GX identified one 6,055-bp scaffold (`scaffold_61`) as fish-derived and
recommended `EXCLUDE`. It contained no final gene annotation. Its removal
reduced the assembly from 167 to 166 scaffolds without changing the
protein-coding annotation.

Final assembly SHA-256:

```text
16455061212d97b94e3f20abcedee353719dc651026adb39744cdb577ca3b2f6
```

## Mitochondrial assembly

```text
/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa
/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.gff3
```

The mitochondrial genome is maintained separately and contains the 12
protein-coding genes expected for cestode mitochondrial genomes.

## Structural annotation

BRAKER3 was run in ETP mode with RNA and related-cestode protein evidence.
Subsequent refinement corrected 3,707 UTRs with TSA evidence, added 187
high-confidence genes, and corrected 4,340 3-prime ends with TagSeq evidence.

A final protein audit removed 62 unsupported, extremely short GeneMark models.
The filtered release contains:

| Feature | Final count |
|---|---:|
| Protein-coding genes | 9,796 |
| Protein-coding transcripts/protein isoforms | 12,849 |
| Exons | 105,095 |
| CDS features | 104,603 |
| 5-prime UTR features | 4,247 |
| 3-prime UTR features | 7,830 |

Filtered longest-isoform proteome BUSCO (`metazoa_odb12`):

```text
C:80.4%[S:79.9%,D:0.4%],F:7.0%,M:12.6%,n:672
```

## ncRNA annotation

The conservative nuclear ncRNA set contains 1,149 features. AGAT creates one
gene parent for each retained ncRNA during the final merge.

| ncRNA type | Count |
|---|---:|
| tRNA | 352 |
| rRNA | 761 |
| snRNA | 32 |
| snoRNA | 2 |
| miRNA | 1 |
| RNase P RNA | 1 |

The high rRNA count reflects retained tandem rDNA arrays. A broader Rfam
candidate set is kept separately and is not part of the conservative release.
No conservative ncRNA feature overlapped a protein-coding CDS in the release
audit.

## Functional annotation

The unfiltered 12,911-protein set was searched against all three resources;
the resulting evidence was then restricted to the 12,849 retained proteins.

| Evidence source | Proteins |
|---|---:|
| Swiss-Prot best hit | 9,172 |
| eggNOG annotation | 8,177 |
| Integrated InterPro entry after filtering | 10,521 |

Conservative product-name assignment used the hierarchy
`high-confidence Swiss-Prot > concise eggNOG description > InterPro domain > hypothetical protein`.
GO/process sentences and generic Reactome propagation were not used as product
names.

| Product category | Proteins |
|---|---:|
| Assigned conservative product | 10,535 (82.0%) |
| Hypothetical protein | 2,314 (18.0%) |
| Missing product | 0 |
| Remaining manual-review records | 0 |

## Final release files

```text
/work/cyu/annotation/ssol-annotation/results/final_filtered_release/
  Ssolidus_nuclear.annotation.final_candidate.gff3
  Ssolidus_nuclear.protein_coding.release_candidate.gtf
  Ssolidus_nuclear.proteins.filtered62.fa
  Ssolidus_nuclear.filtered62.longest_isoform.fa
  functional_integration/Ssolidus_protein_function_master.tsv
  functional_integration/Ssolidus_product_name_assignment.tsv
```

The comprehensive GFF3 contains 10,945 total gene features:

- 9,796 protein-coding genes;
- 1,149 conservative ncRNA genes.

All sequence IDs occur in the final 166-scaffold assembly. The final coordinate
audit contains no invalid, reversed, or out-of-bounds intervals. AGAT parsed all
246,718 feature records successfully.

The versioned shared scientific release is:

```text
/work/shared/Ssolidus_release_v1.0_preNCBI_20260831
```

## Remaining NCBI submission steps

1. Register or confirm the NCBI locus-tag prefix.
2. Assign final submission `locus_tag`, `transcript_id`, and `protein_id` values.
3. Generate the submission package and run `table2asn`/NCBI validator.
4. Archive accessioned files as the NCBI release.

These are submission and accessioning steps; genome assembly, structural
annotation, ncRNA annotation, functional annotation, and scientific release QC
are complete.
