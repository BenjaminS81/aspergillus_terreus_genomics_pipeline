"""Resolve locus order by exact sequence matching without changing source data."""
from pathlib import Path
import json,hashlib
R=Path(r'D:\Illumina\David_Eisele\terreus_genotyping/results/0_4_mafft_iqtree/Sep28')
O=Path(__file__).resolve().parent/'review_reports/2026-11-10'
def fasta(p):
    h=None;a=[]
    with p.open() as f:
        for l in f:
            if l.startswith('>'):
                if h is not None:yield h,''.join(a)
                h=l[1:].strip();a=[]
            else:a.append(l.strip())
    if h is not None:yield h,''.join(a)
target={h.split()[0]:s for h,s in fasta(R/'concatenated.fa')}
ref=target['11_transcripts'];parts=[];bad=[];amb=[]
for i,p in enumerate(sorted((R/'mafft').glob('*_mafft_renamed.fa'))):
    recs=list(fasta(p));seq=next(s for h,s in recs if h.split()[0]=='11_transcripts')
    pos=ref.find(seq)
    if pos<0:
        bad.append({'file':p.name,'issue':'reference sequence not found'});continue
    if ref.find(seq,pos+1)>=0:amb.append(p.name)
    parts.append({'orthogroup':p.name.split('_mafft')[0],'start_1based':pos+1,'end_inclusive':pos+len(seq),'length':len(seq)})
    for h,s in recs:
        key=h.split()[0]
        if target[key][pos:pos+len(seq)]!=s:bad.append({'file':p.name,'sample':key,'issue':'sequence block mismatch'})
    if (i+1)%1000==0:print('Matched ordered blocks',i+1,flush=True)
parts.sort(key=lambda p:p['start_1based'])
continuous=bool(parts) and parts[0]['start_1based']==1 and parts[-1]['end_inclusive']==len(ref) and all(a['end_inclusive']+1==b['start_1based'] for a,b in zip(parts,parts[1:]))
d={'method':'Exact full-locus string location in 11_transcripts; exact block comparison in all 66 rows. No sequence edits.','parts':parts,'ambiguous_reference_matches':amb,'mismatch_count':len(bad),'mismatch_examples':bad[:30],'continuous_complete_coverage':continuous,'samples':len(target),'length':len(ref)}
(O/'concatenation_order_verification.json').write_text(json.dumps(d,indent=2))
print(json.dumps({k:v for k,v in d.items() if k!='parts'},indent=2))
