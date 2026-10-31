#!/usr/bin/env bash
# Environment: genome_assembly. BAM-derived full candidate set, not only the 700 assembled reads.
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
need "$ASM" "$EXTRACT/mt_candidate_reads.fastq.gz" "$SOL_REF" "$SOL"
mkdir -p "$VAL"
[[ ! -e "$VAL/mt.original.paf" && ! -e "$VAL/mt.rotated.paf" ]] || { echo 'Validation outputs exist; inspect them before rerunning.'; exit 1; }
blastn -task blastn -query "$SOL_REF" -subject "$ASM" -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen' -out "$VAL/solidus_mt_vs_scotti_mt.tsv"
blastn -task blastn -query "$SOL" -subject "$ASM" -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen' -out "$VAL/own_solidus_vs_scotti_mt.tsv"
sort -k4,4nr "$VAL/own_solidus_vs_scotti_mt.tsv" > "$VAL/own_solidus_vs_scotti_mt.sorted.tsv"
python "$HELPER" rotate "$ASM" "$VAL/mt.rotated.fa"
minimap2 -x map-hifi -c --secondary=no -t 16 "$VAL/mt.rotated.fa" "$EXTRACT/mt_candidate_reads.fastq.gz" > "$VAL/mt.rotated.paf" 2> "$VAL/mt.rotated.minimap2.log"
python "$HELPER" junction "$VAL/mt.rotated.paf" | tee "$VAL/junction_support.txt"
minimap2 -x map-hifi -c --secondary=no -t 16 "$ASM" "$EXTRACT/mt_candidate_reads.fastq.gz" > "$VAL/mt.original.paf" 2> "$VAL/mt.original.minimap2.log"
python "$HELPER" repeat "$VAL/mt.original.paf" | tee "$VAL/repeat_support.txt"
# Hard-coded repeat window refers ONLY to the unrotated 15437-bp contig_1.
