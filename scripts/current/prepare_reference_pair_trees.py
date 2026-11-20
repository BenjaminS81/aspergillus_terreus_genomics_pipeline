from pathlib import Path
import json,itertools
BASE=Path(__file__).resolve().parent;OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis';P=BASE/'comparative_visualization/2026-11-11/reference_pair_combinations';P.mkdir(exist_ok=True)
def fasta(p):
    d={};k=None
    for line in p.read_text().splitlines():
        if line.startswith('>'):k=line[1:].split()[0];d[k]=''
        else:d[k]+=line.strip()
    return d
data={g:fasta(OLD/'tree_inputs'/(g+'_all_strains.fasta')) for g in ['BenA','CaM','RPB2']};design=[]
for combo in itertools.combinations(data,2):
    name='_'.join(combo);ids=sorted(set.union(*(set(data[g]) for g in combo)));widths=[len(next(iter(data[g].values()))) for g in combo]
    (P/(name+'.fasta')).write_text(''.join('>'+s+'\n'+''.join(data[g].get(s,'?'*w) for g,w in zip(combo,widths))+'\n' for s in ids))
    start=1;parts=[]
    for g,w in zip(combo,widths):parts.append(f'charset {g} = {start}-{start+w-1};');start+=w
    (P/(name+'.nex')).write_text('#nexus\nbegin sets;\n'+'\n'.join(parts)+'\nend;\n');design.append({'name':name,'genes':list(combo),'tips':len(ids),'sites':sum(widths)})
(P/'design.json').write_text(json.dumps(design,indent=2));print(design)
