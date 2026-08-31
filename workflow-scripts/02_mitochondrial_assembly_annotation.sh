#!/usr/bin/env bash
# Mitochondrial ONT read extraction, assembly, polishing and annotation record.
# mt assembly and annotation
set -euo pipefail

###############################################################################
# Step 1. Extract mitochondrial ONT reads
###############################################################################

WORK=/work/cyu/assembly/worm/mitochondria

REF=$WORK/Ssolidus_mt_ref.fa
READS=/work/cyu/assembly/worm/worm_all.top90.q10.fastq.gz
OUT=$WORK/mt_read_reassembly

mkdir -p "$OUT"

###############################################################################
# Map ONT reads to mitochondrial reference
###############################################################################

minimap2 \
    -t 32 \
    -ax map-ont \
    "$REF" \
    "$READS" \
| samtools sort \
    -@ 16 \
    -o "$OUT/ont_vs_mt_ref.bam"

samtools index "$OUT/ont_vs_mt_ref.bam"

###############################################################################
# Mapping statistics
###############################################################################

samtools flagstat \
    "$OUT/ont_vs_mt_ref.bam" \
    > "$OUT/ont_vs_mt_ref.flagstat.txt"

###############################################################################
# Mean coverage
###############################################################################

samtools depth \
-a \
"$OUT/ont_vs_mt_ref.bam" \
| awk '
{
sum+=$3
if($3>max) max=$3
}
END{
print "Mean depth =",sum/NR
print "Max depth =",max
}' \
> "$OUT/mt_depth.txt"

###############################################################################
# Extract mapped reads
###############################################################################

samtools view \
-F 4 \
"$OUT/ont_vs_mt_ref.bam" \
| cut -f1 \
| sort -u \
> "$OUT/mt.readnames.txt"

seqkit grep \
-f "$OUT/mt.readnames.txt" \
"$READS" \
-o "$OUT/mt_reads.fastq.gz"

###############################################################################
# Downsample to ~1000× coverage
###############################################################################

seqkit head \
-n 6000 \
"$OUT/mt_reads.fastq.gz" \
> "$OUT/mt_reads.15Mb.fastq.gz"

#!/usr/bin/env bash
set -euo pipefail

OUT=/work/cyu/assembly/worm/mitochondria/mt_read_reassembly

flye \
--nano-raw "$OUT/mt_reads.15Mb.fastq.gz" \
--genome-size 15k \
--threads 16 \
--out-dir "$OUT/flye_mt"

#!/usr/bin/env bash
set -euo pipefail

cd /work/cyu/assembly/worm/mitochondria/mt_read_reassembly

medaka_consensus \
-i mt_reads.15Mb.fastq.gz \
-d flye_mt/assembly.fasta \
-o medaka_mt \
-t 16 \
-m r1041_e82_400bps_sup_v5.2.0

cat \
flye_mt/assembly_info.txt

makeblastdb \
-in medaka_mt/consensus.fasta \
-dbtype nucl

blastn \
-query ../Ssolidus_mt_ref.fa \
-db medaka_mt/consensus.fasta \
-outfmt "6 qstart qend sstart send length pident" \
| sort -k3,3n

makeblastdb \
-in medaka_mt/consensus.fasta \
-dbtype nucl

blastn \
-query medaka_mt/consensus.fasta \
-db medaka_mt/consensus.fasta \
-outfmt "6 qstart qend sstart send length pident"



conda activate mitos2

cd /work/cyu/assembly/worm/mitochondria/mt_read_reassembly

DBBASE=/work/cyu/assembly/worm/mitochondria/mitos_db

mkdir -p mitos_result

runmitos \
-i medaka_mt/consensus.fasta \
-o mitos_result \
-R "$DBBASE" \
-r refseq89m \
-c 9 \
--best


aragorn \
-mt \
-gc9 \
medaka_mt/consensus.fasta \
> aragorn.medaka.out


cd /work/cyu/assembly/worm/mitochondria/mt_read_reassembly/mitos_result

result.final.gff
