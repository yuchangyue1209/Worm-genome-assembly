# Worm-genome-assembly

Assembly and annotation notes for the *Schistocephalus solidus* genome.

## 项目总结

这个项目记录了 *Schistocephalus solidus* 从原始 ONT 长读长数据到最终基因组和
基因注释的完整流程。工作可以分为五组：

1. **基因组组装与清理**：过滤 ONT reads，组装、polish、去除冗余单倍型，并将
   核基因组和线粒体基因组分开。
2. **重复序列注释**：建立物种特异的重复序列库并对核基因组进行 soft masking。
3. **RNA-seq 证据准备**：将公共 RNA-seq 数据比对到核基因组，生成剪接位点和
   intron hints。
4. **基因结构预测与证据修正**：运行 BRAKER3，然后使用 TSA 和 TagSeq 修正 UTR、
   补充高可信新基因并修正 3' 端。
5. **发布前质量控制**：运行 FCS-GX 污染筛查，审计异常蛋白和 rDNA，注释核
   ncRNA，并检查蛋白编码基因与 ncRNA 是否冲突。
6. **功能注释与发布**：完成 Swiss-Prot 注释；随后整合 eggNOG、InterPro、保守
   product name 和 NCBI 提交字段，最后生成 release v1.0 与校验值。

### 脚本分组和作用

#### 第一组：组装及核/线粒体分离

| 脚本 | 主要作用 | 主要结果 |
|---|---|---|
| `01_ont_nuclear_assembly.sh` | ONT reads 过滤；hifiasm 组装；Racon 和 Medaka polishing；QUAST/BUSCO 检查；purge_dups 去冗余 | 产生经过 polish 和去冗余的候选基因组；这是过程记录，部分命令需要按实际目录分段运行 |
| `02_mitochondrial_assembly_annotation.sh` | 从 ONT reads 中提取线粒体 reads；Flye 重组装；Medaka polish；用 MITOS2 和 ARAGORN 注释 | 线粒体候选序列及线粒体功能注释 |
| `03_document_mt_removal.sh` | 比较移除线粒体前后的 FASTA，找出被删除 scaffold，重新构建不含线粒体的 FASTA并核对序列 | 被移除 scaffold 的 ID/FASTA，以及可重复验证的核基因组 |

#### 第二组：重复序列

| 脚本 | 主要作用 | 主要结果 |
|---|---|---|
| `04_repeat_annotation.sh` | RepeatModeler2 从头构建重复序列库；RepeatMasker 标注并 soft-mask 核基因组 | masked 核基因组、repeat GFF、分类表和重复比例统计 |

最终核基因组有 **68.64%** 重复序列。这个数值描述基因组组成，但不能单独解释早期
Cell Ranger 的低 mapping；后续 nuclear+mt 测试证明，参考中缺少线粒体序列是主要
原因之一。

#### 第三组：RNA-seq 证据

| 脚本 | 主要作用 | 主要结果 |
|---|---|---|
| `05_rnaseq_evidence.sh` | 为 masked 核基因组建立 STAR index；比对 15 个 PRJNA304161 paired-end libraries；保留 `XS:A:` 剪接方向标签；运行 `bam2hints` | 排序并建立索引的 RNA-seq BAM，以及供 BRAKER3 使用的 intron hints |

#### 第四组：BRAKER3 基因预测及初步 QC

| 脚本 | 主要作用 | 主要结果 |
|---|---|---|
| `06_braker3_etp.sh` | 合并 TSA 与近缘物种蛋白证据，过滤异常或过长蛋白，结合 RNA-seq 和蛋白证据运行 BRAKER3 ETP | 原始 BRAKER3 GFF3、GTF、蛋白和 CDS |
| `07_annotation_postprocess_qc.sh` | 用 AGAT 每个基因保留最长蛋白 isoform；检查内部终止密码子；运行 protein-mode BUSCO | BRAKER3 canonical protein set 和中间版本 BUSCO 结果 |

`07_annotation_postprocess_qc.sh` 对应的是 **9,671 genes 的 recovered BRAKER3
中间版本**，不是最终的 9,858-gene release。README 保留该结果，是为了记录最终
注释相对于原始 BRAKER3 的改进基线。

#### 第五组：证据修正和最终发布检查

BRAKER3 之后完成的 TSA/TagSeq 修正产生最终注释：TSA 修正 3,707 个 UTR，加入
187 个高可信新基因，TagSeq 修正 4,340 个 3' ends。对应的最终 GTF 和 GFF3 在
下面的 “Final evidence-refined annotation release” 一节列出。

| 脚本 | 主要作用 | 主要结果 |
|---|---|---|
| `08_structural_release_manifest.sh` | 检查五个最终结构注释文件是否存在；计算 SHA-256；重新统计核基因组长度、scaffold、N50、最大 scaffold 和 GFF3 feature 数量 | 可保存的结构注释 manifest/QC 输出 |

`08_structural_release_manifest.sh` **不会重新生成或覆盖组装和注释**；它只验证已经完成的
final release，适合发布、复制或归档前运行。

### Final release 是什么？

`final_release` 指本项目分析和后续单细胞参考应使用的冻结版本，而不是中间过程文件：

| 文件 | 内容和用途 |
|---|---|
| `/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa` | FCS-GX 清理后的最终核基因组；166 scaffolds，846,268,314 bp；用于核基因组分析和核基因注释 |
| `/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa` | 单独组装的最终线粒体基因组；加入单细胞参考以捕获 mt reads |
| `/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.gff3` | 线粒体基因组的基因注释 |
| `/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.checked.gtf` | 最终核基因 GTF；已经检查为 Cell Ranger 可接受格式 |
| `/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3` | 信息更完整、经 AGAT 规范化的最终核基因 GFF3；用于发布、浏览和常规基因组分析 |

单细胞分析实际使用的是将最终核基因组、线粒体基因组及相应注释合并后建立的
Cell Ranger reference：

```text
/work/cyu/scRNA/reference_test/ssol_final_nuclear_plus_mt
```

简单来说：FASTA 文件回答“基因组序列是什么”，GFF3/GTF 回答“基因在哪里、结构
是什么”；nuclear+mt Cell Ranger reference 则是为 10x reads 定量准备好的索引版本。

## Final nuclear assembly

The release nuclear genome is:

```text
/work/cyu/assembly/worm/final_release/Ssolidus_nuclear.fa
```

The separately assembled mitochondrial genome and its annotation are:

```text
/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.fa
/work/cyu/assembly/worm/final_release/Ssolidus_mitochondrial.gff3
```

### Final assembly metrics

| Metric | Nuclear assembly |
|---|---:|
| Total length | 846,268,314 bp |
| Scaffolds | 166 |
| Scaffold N50 | 99,718,931 bp |
| Longest scaffold | 154,525,006 bp |
| GC content | 44.07% |
| Ambiguous bases (N) | 5,122 bp |
| Repetitive sequence | 68.64% |

FCS-GX identified one 6,055-bp fish-derived scaffold (`scaffold_61`, assigned
to *Esox lucius*) and recommended `EXCLUDE`. It contained no records in the
final nuclear GFF3 or Cell Ranger GTF. Removing it therefore changed the
assembly from 167 to 166 scaffolds without changing the 9,858 protein-coding
genes, proteome, Swiss-Prot results, or completed single-cell quantification.
The current nuclear FASTA SHA-256 is:

```text
16455061212d97b94e3f20abcedee353719dc651026adb39744cdb577ca3b2f6
```

## Workflow scripts

- `01_ont_nuclear_assembly.sh`: ONT filtering, assembly, polishing and purge_dups notes.
- `02_mitochondrial_assembly_annotation.sh`: mitochondrial read extraction, Flye/Medaka
  assembly and MITOS2/ARAGORN annotation.
- `03_document_mt_removal.sh`: identify and archive scaffolds absent from the
  final nuclear FASTA, and reproduce the mt-free assembly.
- `04_repeat_annotation.sh`: RepeatModeler2 and RepeatMasker soft masking.
- `05_rnaseq_evidence.sh`: STAR indexing/alignment and intron hints.
- `06_braker3_etp.sh`: prepare protein evidence and run BRAKER3 with RNA and
  protein evidence (ETP mode).
- `07_annotation_postprocess_qc.sh`: retain the longest protein isoform per
  gene, validate the final proteins and run protein-mode BUSCO.
- `08_structural_release_manifest.sh`: verify the structural release files and report
  checksums, assembly statistics and final GFF3 feature counts.
- `09_fcs_gx_contamination_audit.sh`: document/run FCS-GX and verify the single
  6,055-bp fish contaminant has been removed without deleting annotations.
- `10_final_protein_and_swissprot_audit.sh`: audit unusually short/long final
  proteins and summarize the completed Swiss-Prot DIAMOND annotation.
- `11_swissprot_diamond.sh`: build the Swiss-Prot DIAMOND database, retain the
  top five matches, select one best-bitscore hit per final protein and report
  annotation coverage.
- `12_ncrna_primary_scans.sh`: record the expensive primary tRNAscan-SE,
  barrnap and whole-genome Rfam scans; these are opt-in so completed scans are
  not accidentally repeated.
- `13_ncrna_annotation_and_audit.sh`: reproduce the conservative tRNA/Rfam
  release and test all ncRNA features against protein-coding CDS.
- `14_release_v1_manifest.sh`: validate the current clean assembly, structural
  annotation, ncRNA files and completed Swiss-Prot output, and write checksums.

## Existing annotation outputs on stickleback

```text
/work/cyu/annotation/ssol-annotation/results/repeats
/work/cyu/annotation/ssol-annotation/results/rnaseq
/work/cyu/annotation/ssol-annotation/results/evidence
/work/cyu/annotation/ssol-annotation/results/braker3_recovered
```

The first BRAKER3 run encountered a pathological Spaln comparison between
`Egra_ECG_04546` (4359 aa) and `MSTRG.506` (15190 bp). The updated BRAKER
script excludes that single protein and proteins longer than 10,000 aa when
preparing evidence for a clean rerun.

That run also failed in GeneMark-ETP with `parse_set.pl` phase counts equal to
zero. The original STAR BAM had been made without
`--outSAMstrandField intronMotif`, so spliced alignments lacked the `XS:A:`
strand tags required by StringTie. `05_rnaseq_evidence.sh` now adds the option
and checks that an XS-tagged spliced alignment exists before BRAKER is run.

## Final BRAKER3 annotation

The original run was recovered from its completed GeneMark-ETP results after
patching `getAnnoFastaFromJoingenes.py` for Biopython 1.86 compatibility. The
successful run ended with `BRAKER RUN FINISHED` and produced:

```text
/work/cyu/annotation/ssol-annotation/results/braker3_recovered/braker.gff3
/work/cyu/annotation/ssol-annotation/results/braker3_recovered/braker.gtf
/work/cyu/annotation/ssol-annotation/results/braker3_recovered/braker.aa
/work/cyu/annotation/ssol-annotation/results/braker3_recovered/braker.codingseq
```

Final structural-annotation counts:

```text
genes                         9,671
transcripts/protein isoforms 12,724
proteins with internal stops      0
proteins without terminal stop    8
```

AGAT removed 3,053 shorter isoforms to create a 9,671-protein canonical set:

```text
/work/cyu/annotation/ssol-annotation/results/braker3_recovered/braker.longest_isoform.gff3
/work/cyu/annotation/ssol-annotation/results/braker3_recovered/braker.longest_isoform.aa
```

BUSCO v6.0.0 with `metazoa_odb12` (2026-05-22) recovered
`C:79.6%[S:79.2%,D:0.4%],F:6.1%,M:14.3%,n:672` from the canonical proteins.
The low duplication in the canonical set shows that the 13.8% duplication in
the all-isoform protein run was caused by alternative isoforms. Genome-mode
Miniprot and MetaEuk results were lower and predictor-dependent, so the
canonical protein result is the primary annotation-completeness statistic.

## Final evidence-refined annotation release

The release annotation extends the recovered BRAKER3 models with TSA-supported
UTRs and genes and TagSeq-supported 3' ends. TSA evidence corrected 3,707 UTRs,
187 high-confidence new genes were added, and TagSeq corrected 4,340 3' ends.

```text
# Cell Ranger-compatible nuclear GTF
/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.checked.gtf

# Final AGAT-normalized nuclear GFF3
/work/cyu/annotation/ssol-annotation/results/tsa_alignment/ssol_braker_tsa_utr_tagseq_plus_187.final.agat.gff3
```

| Feature | Count |
|---|---:|
| Genes | 9,858 |
| Transcripts | 12,911 |
| Exons | 105,217 |
| CDS features | 104,725 |
| 5' UTR features | 4,247 |
| 3' UTR features | 7,830 |

Final proteome BUSCO (`metazoa_odb12`):
`C:80.4%[S:79.9%,D:0.4%],F:7.0%,M:12.6%`.

The earlier recovered-BRAKER3 counts and BUSCO result above are retained as an
intermediate baseline; they are not the final release annotation.

## Interpretation relevant to single-cell mapping

Although 68.64% of the nuclear assembly is repetitive, repeats alone did not
cause the low Cell Ranger mapping rate. Adding the separately assembled
mitochondrial genome increased transcriptome mapping by 15.6 percentage points
in Big and 21.1 points in Breeding. Excluded mitochondrial sequence was
therefore a major cause; repeat ambiguity is secondary. A previous minimap2
diagnostic in which 90.97% of low-MAPQ alignments overlapped any repeat must not
be interpreted as evidence of massive LINE transcription or as the primary
Cell Ranger failure mechanism.

## Release-preparation QC completed after structural annotation

### Contamination and excluded rDNA audit

NCBI FCS-GX was run against the final nuclear assembly (`tax-id 70667`, asserted
division `anml:worms`). It found one contaminant:

```text
scaffold_61  1  6055  6055  EXCLUDE  anml:fishes  89  Esox lucius
```

The scaffold had no protein-coding annotation and was removed. Earlier in the
assembly workflow, 26 short rDNA-associated contigs (82,735 bp total) had also
been excluded after depth/rRNA screening. The most substantial excluded contig,
`ptg001004l` (28,068 bp), had 18.31-times median depth and Rfam confirmed a
redundant high-copy rDNA-array fragment (LSU/5.8S/SSU units). It was not restored
because the retained assembly contains complete rDNA arrays, including about 19
units in the first 250 kb of `scaffold_9`.

### Final-protein audit and Swiss-Prot

The final proteome contains 12,911 sequences (7,894,936 aa total; mean 611.5 aa;
maximum 8,390 aa) and no internal stop codons. There are 77 proteins shorter
than 50 aa and 29 longer than 5,000 aa. Of these, one short and 26 long proteins
have Swiss-Prot hits. Sixty-two short models lacking the evidence tags examined
so far remain review candidates; 14 short models have TagSeq support and must
not be removed automatically.

DIAMOND against UniProtKB/Swiss-Prot release 2026_02 produced hits for
9,172/12,911 proteins (71.04%). These results remain valid after FCS cleaning
because the removed scaffold contained no gene models.

### Nuclear ncRNA annotation

The conservative nuclear ncRNA release is:

```text
/work/cyu/annotation/ssol-annotation/results/ncrna_final/release/Ssolidus_nuclear.ncRNA.v1.0.gff3
```

It contains 1,149 features: 761 near-complete rRNA candidates, 352 conservative
tRNAs, 32 snRNAs, 2 snoRNAs, 1 miRNA (`mir-71`) and 1 RNase P RNA. The high rRNA
count reflects tandem high-copy arrays and should not be described as 761
independent functional genes. A broader 2,612-feature Rfam candidate set is
retained separately as supplementary annotation. Neither same-strand nor
opposite-strand intersections were detected between the 1,149 conservative
ncRNA features and protein-coding CDS.

### Current release status

Completed: assembly, polishing, haplotig purging, mt separation/assembly,
repeat annotation, RNA/protein-supported BRAKER3 prediction, TSA/TagSeq
refinement, FCS-GX cleaning, protein audit, Swiss-Prot annotation, rDNA audit
and conservative nuclear ncRNA annotation.

Still required before the publication/NCBI functional release: resolve the 62
short low-confidence protein models, run eggNOG-mapper and InterProScan,
integrate functional evidence, assign conservative product names, add stable
`locus_tag`/`transcript_id`/`protein_id`, run table2asn/NCBI validation, and
freeze the final release with checksums.

The scripts use numeric prefixes to document dependency order. Expensive
primary-scan scripts should not be rerun when their input and database versions
are unchanged; existing completed output should instead be retained and audited.

## Shared annotation files

```text
/work/shared/cyu/ssol-genome-annotation/results/braker3/braker.gff3
/work/shared/cyu/ssol-genome-annotation/results/braker3/braker.aa
/work/shared/cyu/ssol-genome-annotation/results/canonical/braker.longest_isoform.aa
/work/shared/cyu/ssol-genome-annotation/results/evidence/braker_proteins.filtered.fa
```
