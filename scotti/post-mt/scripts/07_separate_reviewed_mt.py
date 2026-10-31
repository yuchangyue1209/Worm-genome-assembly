#!/usr/bin/env python3
"""Reconstruct reviewed RagTag-only mt separation. Never overwrites historical outputs.
Biopython required. Exclusion uses four reviewed primary IDs, not a BLAST threshold.
Usage: python 07_separate_reviewed_mt.py OUTPUT_DIR
"""
from pathlib import Path
from Bio import SeqIO
import sys, hashlib
root=Path('/work/cyu/assembly/worm/scotti')
source=root/'raw/scotti_2_ragtag_500min_direct.fasta'
out=Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=False)
expected={'ptg000539l':30235,'ptg000538l':30810,'ptg000786l':15827,'ptg000644l':12885}
records=list(SeqIO.parse(source,'fasta'))
if (len(records),sum(map(len,records)))!=(416,991406803): raise ValueError('Unexpected assembly')
if len({r.id for r in records})!=len(records): raise ValueError('Duplicate FASTA IDs')
excluded=[r for r in records if r.id in expected]
if {r.id:len(r) for r in excluded}!=expected: raise ValueError('Reviewed mt IDs/lengths differ')
nuclear=[r for r in records if r.id not in expected]
assert (len(nuclear),sum(map(len,nuclear)))==(412,991317046)
SeqIO.write(nuclear,out/'Scotti_nuclear.mt_filtered.ragtag.provisional.fa','fasta')
SeqIO.write(excluded,out/'scotti_mt_candidates.ragtag.fa','fasta')
with open(out/'ragtag_mt_filter_manifest.tsv','x') as f:
    f.write('sequence\tlength_bp\taction\tsequence_sha256\n')
    for r in records:
        action='exclude_reviewed_mt_candidate' if r.id in expected else 'retain_nuclear_candidate'
        f.write(f'{r.id}\t{len(r)}\t{action}\t{hashlib.sha256(str(r.seq).upper().encode()).hexdigest()}\n')
print('Nuclear: 412 sequences; 991,317,046 bp. Excluded: 4; 89,757 bp.')
