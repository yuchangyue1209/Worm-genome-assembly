#!/usr/bin/env bash
# Environment: genome_assembly. Does not use scotti_1 BAM or the malformed old FASTQ.
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
need "$BAM"
TMP="$ROOT/reads/hifi/scotti_2.from_bam.partial.fastq.gz"
[[ ! -e "$READS" && ! -e "$TMP" ]] || { echo 'Output exists; do not overwrite.'; exit 1; }
samtools quickcheck -u -v "$BAM"
samtools view -H "$BAM" > "$ROOT/reads/hifi/scotti_2.bam_header.sam"
samtools fastq -@ 8 -n "$BAM" 2> "$ROOT/reads/hifi/scotti_2.bam_to_fastq.log" | gzip -c > "$TMP"
gzip -t "$TMP"
# Full parse before promoting the partial file; expected count from the observed BAM export.
python "$HELPER" fastq-check "$TMP" --expected 4662246 | tee "$ROOT/reads/hifi/scotti_2.fastq_validation.txt"
mv "$TMP" "$READS"
cat "$ROOT/reads/hifi/scotti_2.bam_to_fastq.log"
