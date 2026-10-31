#!/usr/bin/env bash
# Environment: mitos2 with Biopython, blastn, curl. Evidence only; NO GFF mutation.
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
need "$ASM" "$SOL" "$SOL_GFF"
mkdir -p "$VAL"
python "$HELPER" extract-gff "$SOL" "$SOL_GFF" "$VAL/solidus_rrnL.fa" --genes rrnl
python "$HELPER" extract-gff "$SOL" "$SOL_GFF" "$VAL/solidus_trnH_cox3.fa" --genes trnh cox3
blastn -task blastn -query "$VAL/solidus_rrnL.fa" -subject "$ASM" -evalue 1e-5 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen' -out "$VAL/solidus_rrnL_vs_scotti.tsv"
blastn -task blastn -word_size 7 -dust no -query "$VAL/solidus_trnH_cox3.fa" -subject "$ASM" -evalue 1e-3 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen' -out "$VAL/solidus_trnH_cox3_vs_scotti.tsv"
# Download only if absent; parser below rejects an invalid response or unexpected species.
if [[ ! -s "$VAL/MW602521.1.gb" ]]; then
    curl -fL --retry 3 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MW602521.1&rettype=gbwithparts&retmode=text' -o "$VAL/MW602521.1.gb.partial"
    python "$HELPER" extract-genbank "$VAL/MW602521.1.gb.partial" "$VAL/published_trnH_cox3.fa"
    mv "$VAL/MW602521.1.gb.partial" "$VAL/MW602521.1.gb"
else
    python "$HELPER" extract-genbank "$VAL/MW602521.1.gb" "$VAL/published_trnH_cox3.fa"
fi
blastn -task blastn -word_size 7 -dust no -query "$VAL/published_trnH_cox3.fa" -subject "$ASM" -evalue 1e-3 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen' -out "$VAL/published_trnH_cox3_vs_scotti.tsv"
python "$HELPER" translation "$ASM" "$SOL" | tee "$VAL/cox3_boundary_translation.txt"
cat "$VAL/solidus_rrnL_vs_scotti.tsv"
cat "$VAL/published_trnH_cox3_vs_scotti.tsv"
