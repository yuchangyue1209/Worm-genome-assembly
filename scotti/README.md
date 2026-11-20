# Schistocephalus cotti (scotti_2)

The cotti workflow is independent of [solidus](../solidus/README.md).

- [Mitochondrial workflow](mitochondria/README.md)
- [Post-mitochondrial assembly and annotation scripts](post-mt/README.md)
- **[Current validated working candidate](post-mt/filtered7_working_candidate.md)**
- [Detailed execution history](post-mt/current_status.md)
- [Optional BUSCO model-review investigation](post-mt/candidate_review.md)

Current working candidate: **9,234 protein-coding genes, 11,464 isoforms, 9,234 representative proteins**. Seven reviewed 2–10 aa predictions were excluded without changing genome coordinates or merging neighboring models. GFF3 validation passed; metazoa_odb12 BUSCO remained C73.5%, F9.8%, M16.7%, n672.
Swiss-Prot, eggNOG and InterProScan are complete and integrated for the retained IDs. This is a working candidate, not a final scientific/submission release. TE-related model review, ncRNA annotation, remaining exceptional-model review, and final release verification remain pending.

Selected nuclear genome: RagTag-only, four whole mt candidates removed; 412 sequences, 991,317,046 bp. New SAMBA results remain a separate evaluated experiment. No additional purge_dups removal was performed.
The independent cotti mt working package contains a 15,437 bp assembly and 36 genes. Its revised server v2 curation script is not yet synchronized locally; do not overwrite it with the older local script 08.

Only local code/documentation was organized. No server jobs were submitted by this update.
