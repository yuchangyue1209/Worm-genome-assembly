#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/common.sh"
OUT=${FUNC}/annotation_completeness
mkdir -p "${OUT}"

printf 'Species\tTotal_proteins\tInterPro_any\tInterPro_integrated\tPfam\tInterPro_GO\teggNOG\teggNOG_description\teggNOG_GO\tKEGG_KO\tKEGG_pathway\n' \
    > "${OUT}/functional_coverage.tsv"

for sp in "${SPECIES[@]}"; do
    protein=${STD}/proteins_primary/${sp}.faa
    ipr=${FUNC}/interpro/${sp}.tsv
    egg=${FUNC}/eggnog/${sp}.emapper.annotations
    require_file "${protein}"; require_file "${ipr}"; require_file "${egg}"
    total=$(grep -c '^>' "${protein}")
    read -r any integrated pfam ipr_go < <(awk -F'\t' '
        !a[$1]++{any++}
        $12~/^IPR/&&!i[$1]++{integrated++}
        $4=="Pfam"&&!p[$1]++{pfam++}
        NF>=14&&$14!="-"&&$14!=""&&!g[$1]++{go++}
        END{print any+0,integrated+0,pfam+0,go+0}
    ' "${ipr}")
    read -r egg_all description egg_go ko pathway < <(awk -F'\t' '
        $0!~/^#/ {
            if(!a[$1]++)all++
            if($8!="-"&&$8!=""&&!d[$1]++)description++
            if($10!="-"&&$10!=""&&!g[$1]++)go++
            if($12!="-"&&$12!=""&&!k[$1]++)ko++
            if($13!="-"&&$13!=""&&!p[$1]++)pathway++
        }
        END{print all+0,description+0,go+0,ko+0,pathway+0}
    ' "${egg}")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${sp}" "${total}" "${any}" "${integrated}" "${pfam}" "${ipr_go}" \
        "${egg_all}" "${description}" "${egg_go}" "${ko}" "${pathway}" \
        >> "${OUT}/functional_coverage.tsv"
done

echo "PASS: ${OUT}/functional_coverage.tsv"
