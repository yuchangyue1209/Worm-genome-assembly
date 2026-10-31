# Diagnostic history and unsuccessful approaches

This is a historical record, not a recommended sequence of reruns.

## Malformed original FASTQ

`scotti_2.fastq` contained 52,346,862,294 bytes. SHA-256:
`574fd37221ee67423a9e4018f6d23fa8fc2a96931e47dfb401c656824c9fdb9c`.
The first two records were valid. A bounded record-by-record check encountered
a sequence line exceeding 10 MB at record 4,226, byte offset 75,968,056.
minimap2 2.28-r1209 failed with assertion `ks->seq.l <= (2147483647)`.
The subsequent samtools missing-header error was a consequence of that upstream failure.
The cause was not established as a sequencing, transfer, or source-file generation
error. This does not demonstrate that the original hifiasm input was also malformed.

Historical read-only diagnostics:

```bash
READS=/work/cyu/assembly/worm/scotti/reads/hifi/scotti_2.fastq
stat -c '%n : %s bytes' "$READS"
file "$READS"
sha256sum "$READS"
minimap2 --version
# The helper stops at a line exceeding 10 MB rather than reading indefinitely:
python scripts/mt_helpers.py fastq-check "$READS"
```

Recovery used a new export from `scotti_2_merged.bam`, rather than skipping bad
records. This BAM combines two runs from the same individual; the scotti_1 BAM
was not included.

## Unmapped BAM validation

Plain `samtools quickcheck -v` reported `had no targets in header`.
For unmapped HiFi BAM without a reference dictionary, use
`samtools quickcheck -u -v`.
This checks the header and EOF only. Complete export and full FASTQ parsing
provide additional integrity checks.
BAM may contain kinetics and base-modification information, so the approximately
82G BAM, 52G uncompressed FASTQ, and 6G compressed FASTQ cannot be judged for
corruption by comparing file sizes alone.

## Flye failure with the full candidate set

Historical command:

```bash
# Historical record only; do not use as the default next step.
flye --pacbio-hifi /work/cyu/assembly/worm/scotti/mitochondria/01_read_extraction_from_bam/mt_candidate_reads.fastq.gz --out-dir /work/cyu/assembly/worm/scotti/mitochondria/02_flye_hifi --threads 16
```

On 2026-10-21: total read length 49,450,194 bp; N50/N90 4718/2656;
minimum overlap 3000; overlap-based coverage 2121; zero disjointigs.
A random subset of 700 reads assembled successfully. High coverage was a
troubleshooting hypothesis, not independently established as the sole cause.
The input mode was not changed to ONT, and Medaka was not run.

## Annotation inspection commands

```bash
# Bash; all paths below refer to the server.
OUT=/work/cyu/assembly/worm/scotti/mitochondria/05_mitos2
find "$OUT" -maxdepth 2 -type f -name '*.gff*' -print
cat "$OUT/aragorn.scotti.out"
grep -n -F -B 40 -A 4 'mtRNA-His' "$OUT/aragorn.scotti.out"
# Keep the awk condition on one line to avoid the earlier newline parsing error.
find "$OUT" -maxdepth 2 -type f -name '*.gff*' -exec awk -F '\t' '/^##FASTA/{exit} /^#/{next} NF>=9 && (($4<=5700 && $5>=5000) || tolower($9)~/rrnl|rrns/ || $3=="rRNA" || $3=="origin_of_replication") {print FILENAME "\t" $0}' {} \;
```

The server lacks rg; grep is sufficient. Copy commands from code blocks without
chat-formatting escape characters such as backslashes before asterisks or tildes.
ARAGORN prints each structure before its corresponding name. A structure after
the His label was initially mistaken for His but actually belonged to Gly.
ARAGORN's D-loop tRNA terminology does not identify the mitochondrial control region.

## Evolution of the cox3 boundary interpretation

1. MITOS2 5366–6016: 651 bp, ATG start, 216 amino acids plus TAA; complete ORF but a 63-bp overlap with His.
2. Alignment to the project's solidus annotation, 5372–6016: 645 bp, missing the terminal Q and stop; this range was not adopted.
3. Alignment to the independent published reference, 5428–6070: 643 bp, GTG start, 214 amino acids plus terminal T; consistent with the reference protein span.
4. The third range is the current revision candidate. No final GFF has been generated. Stop completion through polyadenylation is inferred from homology, not demonstrated experimentally in cotti transcripts.

## Outside this directory's scope

FCS-GX reported zero contaminants and zero actions; this does not prove the
absence of all contamination.
Mt-like segments within nuclear scaffolds, the SAMBA ONT-parameter issue,
final nuclear release preparation, and RNA evidence belong to the nuclear
workflow and are not automatically handled here.
