#!/usr/bin/env bash
# Environment: genome_assembly. Use validated BAM-derived FASTQ only.
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
need "$READS" "$CAND"
fresh "$EXTRACT"
minimap2 -ax map-hifi --secondary=no -t 20 "$CAND" "$READS" 2> "$EXTRACT/minimap2.log" | samtools view -@ 4 -b -F 2308 -o "$EXTRACT/mt_candidates.primary.bam" -
samtools fastq -@ 4 -n "$EXTRACT/mt_candidates.primary.bam" 2> "$EXTRACT/samtools_fastq.log" | gzip -c > "$EXTRACT/mt_candidate_reads.fastq.gz"
gzip -t "$EXTRACT/mt_candidate_reads.fastq.gz"
samtools view -c "$EXTRACT/mt_candidates.primary.bam" | tee "$EXTRACT/read_count.txt"
# No MAPQ filter here: four near-duplicate mt candidates make mapping ambiguous.
# Mapped primary reads are retained in full; nuclear mt-like reads can remain.
