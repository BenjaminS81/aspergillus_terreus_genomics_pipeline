from pathlib import Path
import json,collections,csv,numpy as np
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis';T=P/'tree_inputs'
meta=json.loads((P/'tree_reference_metadata.json').read_text());byid={r['reference_id']:r for r in meta}
types=[d for d in meta if d['type_flag'] and d['role']!='historical_reference']
def readfa(path):
    d={};n=None
    for ln in path.read_text().splitlines():
        if ln.startswith('>'):n=ln[1:];d[n]=''
        else:d[n]+=ln.strip().upper()
    return d
result=[];best_by_id=collections.defaultdict(dict)
for l in ['BenA','CaM','RPB2']:
    seqs=readfa(T/(l+'_all_strains.fasta'));ta=[d for d in types if d['reference_id'] in seqs];arr=np.array([list(seqs[d['reference_id']]) for d in ta]);valid=np.isin(arr,list('ACGT'))
    for ident,seq in seqs.items():
        if not ident.startswith(('LOCAL_','HIST2011_')):continue
        q=np.array(list(seq));mask=valid & np.isin(q,list('ACGT'));den=mask.sum(axis=1);num=((arr==q)&mask).sum(axis=1);matches=[]
        for i,d in enumerate(ta):
            if den[i]>=120:matches.append({'reference_id':d['reference_id'],'species':d['analysis_species'],'shared_ACGT_bases':int(den[i]),'identical_bases':int(num[i]),'identity_percent':round(100*num[i]/den[i],5)})
        matches.sort(key=lambda m:(-m['identity_percent'],-m['shared_ACGT_bases']))
        if not matches:continue
        best=matches[0];ties=[m for m in matches if abs(m['identity_percent']-best['identity_percent'])<1e-5]
        row={'id':ident,'locus':l,'best_identity_percent':best['identity_percent'],'best_reference_ids':','.join(m['reference_id'] for m in ties),'best_reference_species':','.join(sorted({m['species'] for m in ties})),'shared_ACGT_bases':best['shared_ACGT_bases'],'ranked_type_matches':matches}
        result.append(row);best_by_id[ident][l]=row
(P/'nearest_type_reference_matches.json').write_text(json.dumps(result,indent=2))
flat=[{k:v for k,v in r.items() if k!='ranked_type_matches'} for r in result]
with (P/'nearest_type_reference_matches.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(flat[0]),delimiter='\t');w.writeheader();w.writerows(flat)
summary=[]
for ident,genes in best_by_id.items():
    sets=[set(d['best_reference_species'].split(',')) for d in genes.values()];common=set.intersection(*sets)
    summary.append({'id':ident,'available_loci':len(genes),'shared_nearest_type_species':','.join(sorted(common)),'status':'concordant_nearest_type_three_loci' if len(common)==1 and len(genes)==3 else 'partial_or_tied_nearest_types' if common else 'nearest_reference_rankings_differ','per_locus':{l:{k:v for k,v in d.items() if k not in ['ranked_type_matches','id']} for l,d in genes.items()}})
(P/'nearest_type_concordance.json').write_text(json.dumps(summary,indent=2))
print('PROJECT CONCORDANCE',collections.Counter((d['status'],d['shared_nearest_type_species']) for d in summary if d['id'].startswith('LOCAL_')))
print('HISTORICAL MINIMUM IDENTITIES',{l:min(d['best_identity_percent'] for d in result if d['locus']==l and d['id'].startswith('HIST')) for l in ['BenA','CaM']})
