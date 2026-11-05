#!/usr/bin/env bash
# Environment: rnaseq_annotation. Completed historically. Retains ptg000291l and secondary alignments.
source "$(dirname "$0")/00_config.sh"
DIR="$ANNO/results/star/02_rna_alignment"
IN="$DIR/Aligned.sortedByCoord.out.bam"
OUT="$DIR/Scotti_nuclear.rnaseq.bam"
TMP="$DIR/Scotti_nuclear.rnaseq.partial.bam"
[[ ! -e "$OUT" && ! -e "$TMP" ]]
samtools view -h "$IN" | awk 'BEGIN{FS=OFS="\t"}
$1=="@SQ" {for(i=2;i<=NF;i++) if($i=="SN:Scotti_mt") next; print; next}
/^@/ {print; next}
$3!="*" && $3!="Scotti_mt" && $7!="Scotti_mt" {print}' \
 | samtools view -@ 8 -b -o "$TMP" -
samtools quickcheck -v "$TMP"
samtools view -H "$TMP" > "$DIR/nuclear.header.sam"
python3 - "$DIR/nuclear.header.sam" <<'EOF'
import sys
names=[]
for line in open(sys.argv[1]):
    f=line.rstrip().split('\t')
    if f[0]=='@SQ': names.append(next(x[3:] for x in f if x.startswith('SN:')))
assert len(names)==len(set(names))==412
assert 'Scotti_mt' not in names and 'ptg000291l' in names
print('PASS: 412 nuclear references; mt excluded; ptg000291l retained')
EOF
mv "$TMP" "$OUT"
samtools index -@ 8 "$OUT"
samtools flagstat -@ 8 "$OUT" | tee "$DIR/nuclear.flagstat.txt"
# Full scan avoids early-exit SIGPIPE under pipefail.
samtools view "$OUT" | awk -F '\t' '$6~/N/ {n++; for(i=12;i<=NF;i++) if($i~/^XS:A:[+-]$/){xs++;break}} END {print "Spliced records:",n+0; print "With XS:",xs+0}' \
 | tee "$DIR/nuclear.XS_audit.txt"
