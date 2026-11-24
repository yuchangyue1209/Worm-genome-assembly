"""Check reverse-strand conversion, score headers, tags, CDS unions and output."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess,sys,json,csv
here=Path(__file__).resolve().parent
with TemporaryDirectory() as temp:
    root=Path(temp);base=root/'results/ncrna_v1'
    source=base/'03b_trna_highconfidence_breakdown';source.mkdir(parents=True)
    rfam=base/'02_rfam';rfam.mkdir()
    pc=root/'results/filtered_candidates/Scotti_filtered7_v1';pc.mkdir(parents=True)
    genome=root/'genome.fa';genome.write_text('>chr1\n'+'A'*2000+'\n>chr2\n'+'A'*2000+'\n')
    h1=['Sequence','','tRNA','Bounds','tRNA','Anti','Intron Bounds','Inf','HMM',"2'Str",'Isotype','Isotype','']
    h2=['Name','tRNA#','Begin','End','Type','Codon','Begin','End','Score','Score','Score','CM','Score','Note']
    rows=[['chr1',1,101,180,'Ala','AGC',0,0,80,50,30,'Ala',90,'high confidence set'],
          ['chr1',2,1080,1000,'Arg','ACG',0,0,85,55,30,'Arg',95,'high confidence set'],
          ['chr2',1,10,100,'SeC','TCA',0,0,75,50,25,'SeC',80,'high confidence set'],
          ['chr2',2,201,280,'Ala','AGC',0,0,80,50,30,'Arg',90,'IPD:-2,isotype mismatch'],
          ['chr2',3,501,580,'Undet','NNN',0,0,25,20,5,'Leu',30,'pseudo']]
    source.joinpath('Scotti_tRNA.out').write_text('\t'.join(h1)+'\n'+'\t'.join(h2)+'\n'+'--------\n'+''.join('\t'.join(map(str,r))+'\n' for r in rows))
    source.joinpath('Scotti_tRNA.log').write_text('High confidence set: 3\n')
    gff='##gff-version 3\n'
    for start,end,strand in [(100,130,'+'),(120,140,'+'),(160,170,'-'),(1050,1060,'-')]:
        gff+=f'chr1\tfixture\tCDS\t{start}\t{end}\t.\t{strand}\t0\tParent=tx\n'
    pc.joinpath('Scotti_nuclear.filtered7.candidate.gff3').write_text(gff)
    fields=['1','U1','RF00003','chr1','-','CL00001','cm','1','100','16','100','+','no','1','0.5','0.0','100','1e-20','!','^','-','-','-','-','-','-','100','2000','U1 spliceosomal RNA']
    rfam.joinpath('Scotti_nuclear.rfam.tblout').write_text('# format 2\n'+' '.join(fields)+'\n')
    run=subprocess.run([sys.executable,str(here/'build_scotti_trna_candidate.py'),'--root',str(root),'--genome',str(genome),'--expected','3'],capture_output=True,text=True,check=True)
    out=base/'04_trna_candidate_v1';summary=json.loads(out.joinpath('summary.json').read_text())
    assert summary['high_confidence_tRNA_features']==3
    assert summary['isotype_mismatch_review']==1
    assert summary['tRNA_features_overlapping_any_CDS']==2
    assert summary['tRNA_features_overlapping_same_strand_CDS']==2
    with out.joinpath('tRNA_high_confidence.audit.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
    assert int(r[0]['CDS_overlap_bp'])==51 and int(r[0]['same_strand_CDS_overlap_bp'])==40
    assert r[1]['start']=='1000' and r[1]['end']=='1080' and r[1]['strand']=='-'
    features=[l.split('\t') for l in out.joinpath('Scotti_nuclear.tRNA.high_confidence.candidate.gff3').read_text().splitlines() if not l.startswith('#')]
    assert len(features)==3 and all(len(l)==9 for l in features)
    assert features[1][3:5]==['1000','1080']
    bed=out.joinpath('Scotti_nuclear.tRNA.high_confidence.candidate.bed').read_text().splitlines()[1].split('\t')
    assert bed[1:3]==['999','1080']
    assert summary['Rfam_raw_hits']==1
    assert summary['isotypes']=={'Ala':1,'Arg':1,'SeC':1}
print('PASS: score-column parsing; high-confidence selection; reverse-strand GFF/BED; union CDS coverage; mismatch export; Rfam fmt2/genome checks.')
