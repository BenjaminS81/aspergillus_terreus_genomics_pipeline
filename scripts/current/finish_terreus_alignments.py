from pathlib import Path
import json,collections,csv
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis';OUT=P/'tree_inputs';OUT.mkdir(exist_ok=True)
refs=json.loads((P/'tree_reference_metadata.json').read_text());refids={d['reference_id'] for d in refs}
def readfa(p):
    d={};n=None
    for line in p.read_text().splitlines():
        if line.startswith('>'):n=line[1:].split()[0].removeprefix('_R_');assert n not in d;d[n]=''
        else:d[n]+=line.strip().upper()
    return d
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
summary=[];trimmed={};cols={}
for l in ['BenA','CaM','RPB2']:
    path=P/(l+'_with_historical_aligned.fasta') if l in ['BenA','CaM'] else P/(l+'_genomic_combined_aligned.fasta')
    fa=readfa(path);lengths={len(s) for s in fa.values()};assert len(lengths)==1;length=lengths.pop();rf=[s for k,s in fa.items() if k in refids]
    keep=[i for i in range(length) if sum(s[i] in 'ACGT' for s in rf)/len(rf)>=.5]
    clean={k:''.join(s[i] for i in keep) for k,s in fa.items()};trimmed[l]=clean;cols[l]=[i+1 for i in keep]
    writefa(OUT/(l+'_all_strains.fasta'),clean)
    summary.append({'locus':l,'sequences':len(clean),'reference_sequences':len(rf),'project_sequences':sum(k.startswith('LOCAL_') for k in fa),'untrimmed_columns':length,'retained_columns':len(keep),'criterion':'ACGT occupancy >=50% among reference sequences; no imputation'})
ids=sorted(set.union(*(set(d) for d in trimmed.values())));concat={k:''.join(trimmed[l].get(k,'?'*len(cols[l])) for l in trimmed) for k in ids};writefa(OUT/'combined_all_strains.fasta',concat)
focus={d['reference_id'] for d in refs if d['type_flag'] or d['accepted_type_name'] or d['historical_type_name']};focus|={k for k in ids if k.startswith('LOCAL_')}
writefa(OUT/'combined_type_focused.fasta',{k:s for k,s in concat.items() if k in focus})
start=1;parts=[]
for l in trimmed:
    end=start+len(cols[l])-1;parts.append(f'charset {l} = {start}-{end};');start=end+1
(OUT/'partitions.nex').write_text('#nexus\nbegin sets;\n'+'\n'.join(parts)+'\nend;\n')
(OUT/'alignment_summary.json').write_text(json.dumps(summary,indent=2));(OUT/'retained_column_map.json').write_text(json.dumps(cols,indent=2))
occupancy=[]
for k in ids:
    d={'id':k,'role':'project' if k.startswith('LOCAL_') else 'reference'}
    for l in trimmed:d[l+'_ACGT']=sum(c in 'ACGT' for c in trimmed[l].get(k,''))
    d['available_loci']=sum(d[l+'_ACGT']>=120 for l in trimmed);occupancy.append(d)
(OUT/'locus_occupancy.json').write_text(json.dumps(occupancy,indent=2))
print('ALIGNMENTS',summary,'CONCATENATED',len(concat),'TYPE_FOCUSED',len(focus))
