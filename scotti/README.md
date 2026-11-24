# Schistocephalus cotti (scotti_2)

**Scientific release v1.0 (pre-NCBI) was published on 20 November 2026.**

```text
/work/shared/Scotti_release_v1.0_preNCBI_20261118
```

The selected baseline contains **9,234 protein-coding genes, 11,464 protein isoforms and 9,234 longest representative proteins**. All 32 strong TE-associated candidate genes and the 47 unresolved TE cases are retained. The 9,202-protein sensitivity set is provenance only.

- **[Release files, counts and validation](release_v1.0.md)**
- [Mitochondrial workflow and archived candidate annotation](mitochondria/README.md)
- [Post-mitochondrial workflow scripts](post-mt/README.md)
- [ncRNA, snapshot and shared-publication record](post-mt/ncrna_and_release.md)
- [Historical filtered7 source package](post-mt/filtered7_working_candidate.md)
- [Detailed execution history](post-mt/current_status.md)
- [Optional BUSCO model-review investigation](post-mt/candidate_review.md)

Selected nuclear genome: RagTag-only, with four whole mt candidate sequences removed; 412 sequences and 991,317,046 bp. New SAMBA outputs remain a separate evaluated experiment. No additional purge_dups removal was performed.

Swiss-Prot, eggNOG and InterProScan are complete and integrated for the retained IDs. The shared package includes a separate 3,681-locus nuclear ncRNA candidate component and a 15,437-bp mitochondrial genome with 36 genes. Protein-coding, nuclear ncRNA and mitochondrial GFF3 validation passed. Functional evidence is not a curated product-name release, and ncRNA functional/pseudogene status is not universally resolved.

The authoritative server mitochondrial v2 builder is preserved in the shared mitochondrial package. The older local script 08 is historical and must not overwrite it. See the release record for exact paths.

This workflow is independent of [solidus](../solidus/README.md). Local code/documentation updates do not launch server jobs or replace frozen comparative inputs.
