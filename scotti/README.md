# Schistocephalus cotti (scotti_2)

The cotti workflow is maintained separately from [solidus](../solidus/README.md).

- [Mitochondrial assembly and early annotation workflow](mitochondria/README.md)
- **[After mitochondrial assembly: nuclear assembly assessment, repeats, and RNA processing](post-mt/README.md)**
- [Results, assembly selection, and documentation gaps](post-mt/execution_record.md)
- [Historical issues](mitochondria/history_and_caveats.md)

The selected nuclear assembly is the RagTag-only assembly with four whole mitochondrial candidate sequences removed: 412 sequences, 991,317,046 bp.
The new SAMBA assembly was evaluated but was not selected for annotation. No additional purge_dups removal was performed.

The mitochondrial working package on the server contains a 15,437 bp assembly and 36 annotated genes. Its GFF3 passed GenomeTools validation; it remains a working candidate. The revised server-side v2 curation script has not been synchronized locally. Do not overwrite the validated working package using the earlier local script 08.

RepeatModeler, RepeatMasker, full rRNA screening, STAR mapping, and nuclear BAM preparation are complete. BRAKER3 ETP completed successfully: 9,241 genes and 11,471 protein isoforms. Protein QC, BUSCO, and Swiss-Prot are complete. The eggNOG command has been corrected to -i/-m; successful execution has not yet been reported. InterProScan remains pending. See [current status and results](post-mt/current_status.md).

This update organizes local code and documentation only. It does not submit or rerun server analyses.
