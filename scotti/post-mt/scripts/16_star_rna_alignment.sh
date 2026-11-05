#!/usr/bin/env bash
# Environment: rnaseq_annotation. Completed historically; retains nuclear and mt alignments.
source "$(dirname "$0")/00_config.sh"
QC="$ANNO/results/rnaseq_qc/05_sortmerna_full"
INDEX="$ANNO/results/star/01_index_no_gtf"
OUT="$ANNO/results/star/02_rna_alignment"
mkdir "$OUT"
python3 - "$QC" "$OUT" <<'EOF'
from pathlib import Path
import sys
qc,out=map(Path,sys.argv[1:]); lists=[]
for mate in (1,2):
    paths=(qc/f'non_rrna_R{mate}.list').read_text().splitlines()
    if not paths: raise ValueError('Empty input list')
    for p in paths:
        if not Path(p).is_file() or not Path(p).stat().st_size: raise ValueError(p)
    lists.append(paths)
if len(lists[0])!=len(lists[1]): raise ValueError('Mate file counts differ')
for r1,r2 in zip(*lists):
    if Path(r1).name.replace('_R1.','_R2.')!=Path(r2).name: raise ValueError('Mate mismatch')
for mate,paths in enumerate(lists,1): (out/f'R{mate}.inputs.txt').write_text(','.join(paths))
EOF
READ1=$(cat "$OUT/R1.inputs.txt"); READ2=$(cat "$OUT/R2.inputs.txt")
STAR --runThreadN 20 --genomeDir "$INDEX" --readFilesIn "$READ1" "$READ2" \
  --readFilesCommand zcat --twopassMode Basic --sjdbOverhang 100 \
  --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 40000000000 --outFileNamePrefix "$OUT/" 2>&1 | tee "$OUT/run.log"
samtools index -@ 8 "$OUT/Aligned.sortedByCoord.out.bam"
samtools quickcheck -v "$OUT/Aligned.sortedByCoord.out.bam"
cat "$OUT/Log.final.out"
