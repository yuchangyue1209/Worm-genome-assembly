#!/usr/bin/env bash
# Prepare protein evidence and run BRAKER3 in ETP mode using RNA-seq and protein
# evidence. Run in the braker3_annotation Conda environment.
set -euo pipefail

ROOT=/work/cyu/annotation/ssol-annotation
THREADS=32
GENOME=${ROOT}/results/repeats/worm.q10_r500k_no_contig_ec.nuclear.final.fa.masked
BAM=${ROOT}/results/rnaseq/star_align/Aligned.sortedByCoord.out.bam
EVIDENCE=${ROOT}/results/evidence
RAW=${EVIDENCE}/braker_proteins.raw.fa
FILTERED=${EVIDENCE}/braker_proteins.filtered.fa
BRAKER_DIR=${ROOT}/results/braker3

TSA=${ROOT}/../transcriptomes/hebert-2016/PRJNA304161/GEEE01/GEEE01.1.fsa_aa.gz
REFS=(
    "${ROOT}/../genomes/downloads/spro/ncbi_dataset/data/GCA_053814205.1/protein.faa"
    "${ROOT}/../genomes/downloads/dlat/wormbase_parasite/dibothriocephalus_latus.PRJEB1206.WBPS19.protein.fa"
    "${ROOT}/../genomes/downloads/seri/wormbase_parasite/spirometra_erinaceieuropaei.PRJEB1202.WBPS19.protein.fa"
    "${ROOT}/../genomes/downloads/hmic/wormbase_parasite/hymenolepis_microstoma.PRJEB124.WBPS19.protein.fa"
    "${ROOT}/../genomes/downloads/egra/wormbase_parasite/echinococcus_granulosus.PRJNA754835.WBPS19.protein.fa"
)

export AUGUSTUS_CONFIG_PATH=${ROOT}/config/augustus
export GENEMARK_PATH=/work/cyu/annotation/software/GeneMark-ETP/bin
export PROTHINT_PATH=/work/cyu/annotation/software/GeneMark-ETP/bin/gmes/ProtHint/bin
export PATH="${GENEMARK_PATH}:${PROTHINT_PATH}:/work/cyu/annotation/software/GeneMark-ETP/tools:${PATH}"

for cmd in braker.pl augustus gmetp.pl prothint.py seqkit samtools; do
    command -v "${cmd}" >/dev/null || { echo "ERROR: ${cmd} is not in PATH" >&2; exit 1; }
done
for file in "${GENOME}" "${BAM}" "${TSA}" "${REFS[@]}"; do
    [[ -s "${file}" ]] || { echo "ERROR: missing ${file}" >&2; exit 1; }
done
[[ -w "${AUGUSTUS_CONFIG_PATH}/species" ]] || {
    echo "ERROR: AUGUSTUS config species directory is not writable" >&2; exit 1;
}

mkdir -p "${EVIDENCE}" "${ROOT}/logs"
: > "${RAW}"

append_fasta() {
    local label=$1 input=$2
    if [[ "${input}" == *.gz ]]; then gzip -cd "${input}"; else command cat "${input}"; fi |
    awk -v prefix="${label}" '
        /^>/ {sub(/^>/,""); split($0,a,/[ \t]+/); id=a[1]; gsub(/[^A-Za-z0-9_.:-]/,"_",id); print ">" prefix "_" id; next}
        {gsub(/[[:space:]-]/,""); if(length($0)) print toupper($0)}
    ' >> "${RAW}"
}

append_fasta Ssol_TSA "${TSA}"
labels=(Spro Dlat Seri Hmic Egra)
for i in "${!REFS[@]}"; do append_fasta "${labels[$i]}" "${REFS[$i]}"; done

# The original run spent >9 days on Egra_ECG_04546 versus MSTRG.506.
seqkit grep -v -r -p '^Egra_ECG_04546$' "${RAW}" |
    seqkit seq -m 30 -M 10000 -g |
    seqkit rmdup -s -i -j "${THREADS}" -o "${FILTERED}"
seqkit stats "${RAW}" "${FILTERED}"

if [[ -d "${BRAKER_DIR}" ]] && find "${BRAKER_DIR}" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: existing BRAKER directory is not empty: ${BRAKER_DIR}" >&2
    echo "Rename and preserve it before starting a clean rerun." >&2
    exit 1
fi
mkdir -p "${BRAKER_DIR}"

braker.pl --genome="${GENOME}" --bam="${BAM}" --prot_seq="${FILTERED}" \
    --softmasking --species=ssol_cyu_2026 --workingdir="${BRAKER_DIR}" \
    --threads="${THREADS}" --GENEMARK_PATH="${GENEMARK_PATH}" \
    --PROTHINT_PATH="${PROTHINT_PATH}" --gff3
