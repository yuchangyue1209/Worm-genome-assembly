# Results and workflow decisions

## Assembly comparison

|Assembly|Sequences|Bases|N bases|
|---|---:|---:|---:|
|hifiasm primary|972|991,351,203|0|
|RagTag before SAMBA|416|991,406,803|55,600|
|Old SAMBA, ont preset|413|991,407,400|55,600|
|New SAMBA, HiFi asm preset|414|991,408,080|55,600|

New SAMBA added a net 1,277 bp through two merges:

- scaffold_9_RagTag + scaffold_16_RagTag → 39,789,088 bp.
- ptg000897l + ptg000751l → 33,733 bp.

Changed sequences did not match the newly assembled mitochondrial genome. These changes therefore do not show SAMBA joining the previously identified duplicated mitochondrial candidates.

Original-to-new PAF coordinates (zero-based, half-open):

- scaffold_9: query 0–32,128,717 → target 0–32,128,717.
- scaffold_16: query 3–7,659,419 → target 32,129,672–39,789,088.
- ptg000897l: query 0–16,132 → target 0–16,132.
- ptg000751l: query 7–17,255 → target 16,485–33,733.

Full HiFi mapping to new SAMBA: 4,662,246 primary reads; 4,657,066 mapped primary reads (99.89%); 140,930 supplementary alignments.

Local junction results were supplied by the user. The complete original checking script was not retained and has not been reconstructed here:

|Junction|Spanning reads, MAPQ ≥20|Passing all filters|
|---|---:|---:|
|Large scaffold merge|2|1|
|Short contig merge|2|2|

Reported filters included at least 500 bp extending to each side, identity ≥95%, at least 450 aligned bases in each flank, and no indel ≥50 bp.
Sparse support does not establish a misassembly. The conservative current choice is RagTag-only, without adopting the new SAMBA joins.

## Primary assembly purge assessment

The main coverage peak was 26–28×, with 27× the largest bin.

|Depth|Assembly fraction|
|---|---:|
|0×|0.01%|
|1–5×|0.58%|
|6–10×|3.74%|
|11–16×|10.08%|
|17–20×|10.74%|
|21–35×|62.72%|
|36–50×|10.94%|
|51–500×|1.18%|

BUSCO metazoa_odb10: C73.1% [S71.2%, D1.9%], F3.2%, M23.7%, n954, E31.0%. Of the complete BUSCOs, 216 contained internal stops.
No additional purge_dups removal was performed based on these results. Low duplication does not exclude every haplotig, and low coverage alone does not establish redundancy. The internal-stop issue remains unresolved.

## Mitochondrial separation

The same four whole mitochondrial candidate sequences remained across primary, RagTag, and SAMBA:

|ID|Length, bp|
|---|---:|
|ptg000539l|30,235|
|ptg000538l|30,810|
|ptg000786l|15,827|
|ptg000644l|12,885|

These four sequences were removed from the selected RagTag assembly. Local mitochondrial matches in large nuclear contigs already existed in the primary assembly and were retained. Their individual status as NUMTs has not been established.

Keep the following nuclear versions distinct:

- Old SAMBA: 409 sequences, 991,317,643 bp; original provisional file.
- New SAMBA: 410 sequences, 991,318,323 bp; samba_hifi.provisional file.
- **Selected RagTag: 412 sequences, 991,317,046 bp; ragtag.provisional file.**

An earlier FCS-GX run reported zero contaminant actions. This does not establish that every subsequent assembly version was independently screened. FCS-GX was not rerun during this organization task.

## Missing original code and synchronization gaps

The conversation did not retain complete original code for junction support counts, primary depth plotting, mitochondrial interval-union summaries, overrepresented RNA extraction, or construction of the rRNA reference. Their results are recorded from user-supplied output.
The scripts provide reviewable reconstructions of key steps without claiming recovery of missing original commands, software versions, Hi-C evidence, or new server results.

The revised server-side mitochondrial v2 curation script has not been synchronized locally. The local `mitochondria/scripts/08` remains an earlier candidate version; do not use it to overwrite the validated server working package.
