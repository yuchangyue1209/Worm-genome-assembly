#!/usr/bin/env python3
"""Record evidence-based provisional review of the first 30 primary candidates.
Source sequences and the initial screening outputs are never modified.
"""
from pathlib import Path
from collections import Counter
import argparse, csv, hashlib, json, shutil

# All calls are provisional and refer only to the supplied model/evidence.
DECISIONS={}
def add(genes,category,action,reason):
    for gene in genes.split():
        assert gene not in DECISIONS
        DECISIONS[gene]=(category,action,reason)

add('g1102 g7873','strong_TE_associated_candidate','candidate_for_TE_exclusion_comparison_not_executed',
    'Near-complete classified-TE CDS overlap plus Pol homology and RT/integrase/RNase-H architecture. TE-associated model is strongly supported; this does not prove an active element or justify deleting genome sequence.')
add('g4531','strong_TE_associated_candidate','candidate_for_TE_exclusion_comparison_not_executed',
    '100% classified-TE CDS overlap; high query coverage Pol hit and integrase core/zinc-binding evidence. Short TE-associated fragment, not proof of an active full-length element.')
add('g6406','strong_TE_associated_candidate','candidate_for_TE_exclusion_comparison_not_executed',
    '100% classified-TE CDS overlap; high query coverage Pol hit and RT-domain evidence. Very short fragment; keep source model and locus available for audit.')
add('g6446','strong_TE_associated_candidate','candidate_for_TE_exclusion_comparison_not_executed',
    '100% classified-TE CDS overlap; 88.5% query coverage LINE ORF2 hit and RT description. Integrated domain evidence is absent; candidate classification remains provisional.')
add('g520','strong_TE_associated_candidate','candidate_for_TE_exclusion_comparison_not_executed',
    'Tc1-like transposase DDE annotation plus classified-TE/unknown repeat overlap. Very short fragment; do not equate this with an intact autonomous transposon.')

CELLULAR={
'g112':'SET/Pre-SET/Post-SET architecture and high-coverage histone methyltransferase homology; SETMAR best-hit name does not establish a worm SETMAR ortholog or a transposon.',
'g1226':'RECQ-like helicase homology and helicase-domain architecture; a shared zinc-finger signal does not establish TE origin.',
'g1365':'SRSF7-like homology and splicing-factor/RRM annotation; shared CCHC motifs do not establish TE origin.',
'g2456':'Splicing factor 1 homology with KH and SF1-specific annotations; retain as a cellular-function-supported model.',
'g301':'Solute-carrier homology with SPX/EXS domains; no TE-specific architecture is shown in the supplied summary.',
'g3824':'RBBP6-family and DWNN-domain evidence supports a cellular protein; shared zinc-finger annotations alone are insufficient for TE exclusion.',
'g5601':'TRAMP/AT-rich-DNA-binding annotation and CCHC domains; no TE-specific architecture is shown in the supplied summary.',
'g5806':'CPSF4-like homology and family annotation supports RNA processing; shared zinc fingers do not establish TE origin.',
'g6022':'ZCRB1-family/RRM evidence and high-coverage homology supports a cellular splicing-related model.',
'g6936':'DDI1-type peptidase annotation and DDI1 homology; retroviral-like protease fold is shared with cellular DDI proteins and is not sufficient for TE exclusion.',
}
for gene,reason in CELLULAR.items():
    add(gene,'cellular_function_supported','retain_in_current_proteome_pending_full_audit',reason)
add('g4139 g6181 g6980','generic_CCHC_review','retain_pending_exact_signature_review',
    'CCHC/zinc-knuckle evidence is shared by many proteins. The supplied summary lacks TE-specific evidence sufficient for exclusion. Inspect the exact original screening trigger.')
add('g1071 g5909','TE_association_plausible','retain_pending_domain_and_repeat_family_review',
    'High classified-TE CDS overlap plus eggNOG retroelement description, but integrated InterPro and Swiss-Prot evidence is absent. Inspect member-signature evidence and the matching repeat family.')
add('g200 g5344','TE_association_plausible','retain_pending_domain_and_repeat_family_review',
    'eggNOG Pao-retrotransposon peptidase description; classified-TE overlap is partial and integrated domain/Swiss-Prot evidence is absent. Not enough for exclusion yet.')
add('g3525 g4404 g8040','exact_signature_needed','retain_pending_exact_signature_review',
    'Original priority flag is not explained by the displayed Swiss-Prot/eggNOG/integrated-InterPro summary. Inspect the original member-database description, score and coordinates; repeat overlap alone is not sufficient.')
add('g3174','conflicting_RT_annotation','retain_pending_homology_and_domain_review',
    'Telomerase best hit has only 31.8% query coverage and borderline screening E-value, while InterPro includes group-II-intron/mobility annotation and eggNOG gives a discordant kinase description. Neither TERT identity nor TE origin is established.')
add('g4424','cellular_TE_family_conflict','retain_pending_homology_and_domain_review',
    'eggNOG MULE-transposase description conflicts with ZSWIM1/3 family/domain context and zero repeat overlap. Check domain identity and coverage; do not delete on the MULE label alone.')
add('g8791','weak_retroelement_annotation','retain_pending_homology_and_domain_review',
    'eggNOG retroelement description without repeat overlap or a TE-specific integrated domain. DUF7041 alone does not resolve TE origin.')
add('g9033','possible_mixed_model','retain_pending_domain_layout_and_model_review',
    '1190-aa model combines MIZ/SP-RING and transposase-family annotations; SUMO-ligase Swiss-Prot hit covers only 7.1%. Inspect domain positions and repeat/CDS overlap for a fusion, fragmented model or other explanation; do not label the whole protein from the short best hit.')
assert len(DECISIONS)==30

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--review',type=Path,default=Path('/work/cyu/annotation/scotti-annotation/results/annotation_qc/08_TE_protein_review_v1'))
    p.add_argument('--expected-primary',type=int,default=9234)
    a=p.parse_args();src=a.review
    out=src/'first30_case_review_v1'
    if out.exists():raise FileExistsError(f'Output exists: {out}')
    source=src/'Scotti_TE_review.primary.tsv'
    with source.open() as f:
        reader=csv.DictReader(f,delimiter='\t');fields=list(reader.fieldnames);rows=list(reader)
    assert len(rows)==a.expected_primary
    assert len({r['gene_id'] for r in rows})==len(rows)
    selected=[r for r in rows if r['gene_id'] in DECISIONS]
    assert {r['gene_id'] for r in selected}==set(DECISIONS)
    assert all(r['review_tier'] in {'priority_functional_review','host_context_review'} for r in selected)
    assert sum(r['review_tier']=='host_context_review' for r in selected)==1
    extra=['case_review_category','provisional_action','case_review_reason']
    for r in selected:
        category,action,reason=DECISIONS[r['gene_id']]
        r.update(zip(extra,(category,action,reason)))
    out.mkdir()
    with (out/'first30.primary_case_review.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t');w.writeheader();w.writerows(selected)
    strong=[r for r in selected if r['case_review_category']=='strong_TE_associated_candidate']
    (out/'strong_TE_associated.primary_candidate_ids.txt').write_text(''.join(r['protein_id']+'\n' for r in strong))
    # Preserve exact triggering descriptions and domain coordinates for the next review.
    selected_ids={r['protein_id'] for r in selected}
    match_source=src/'matched_interpro_evidence.tsv'
    with match_source.open() as f, (out/'first30.matched_interpro_evidence.tsv').open('x') as dst:
        reader=csv.DictReader(f,delimiter='\t');w=csv.DictWriter(dst,fieldnames=reader.fieldnames,delimiter='\t');w.writeheader()
        for r in reader:
            if r['protein_id'] in selected_ids:w.writerow(r)
    compact_fields=['gene_id','protein_id','length_aa','CDS_RM_classified_TE_pct','explicit_TE_signals','RT_integrase_signals','generic_signals','host_context_signals','case_review_category']
    with (out/'first30.screening_triggers.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=compact_fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(selected)
    remaining=[r for r in rows if r['review_tier']=='functional_candidate_review']
    with (out/'remaining_functional_candidates.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(remaining)
    summary={'status':'provisional_case_review_no_models_removed','reviewed_primary_models':len(selected),
        'category_counts':dict(Counter(r['case_review_category'] for r in selected)),
        'strong_TE_associated_candidates':len(strong),'remaining_functional_candidates':len(remaining)}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'README.txt').write_text('First 30 primary candidates, based on user-supplied evidence. Provisional actions; no sequences or annotation features removed. Six strong TE-associated candidates are retained in the source proteome and listed only for a potential exclusion sensitivity comparison, not as a validated final deletion set. Specific InterPro member descriptions/coordinates and original triggers are exported for checking. Ten cellular-function-supported models are provisionally retained; exact orthology/product identity is not assigned from best-hit names. Generic CCHC and retroviral-like peptidase descriptors are not sufficient TE evidence. 70 functional candidates remain unreviewed if input is the reported complete screening table. This is a case-specific overlay, not a rewrite of global screening rules.\n')
    shutil.copy2(Path(__file__),out/'review_scotti_te_first30.py')
    with (out/'input.SHA256SUMS').open('x') as f:
        for path in [source,match_source]:f.write(f'{sha(path)}  {path}\n')
    with (out/'SHA256SUMS').open('x') as f:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.name!='SHA256SUMS':f.write(f'{sha(path)}  {path.name}\n')
    print(json.dumps(summary,indent=2));print('Output:',out)

if __name__=='__main__':main()
