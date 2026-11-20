from pathlib import Path
import json,itertools
BASE=Path(__file__).resolve().parent;OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis';P=BASE/'comparative_visualization/2026-11-11/core66_combinations';P.mkdir(exist_ok=True)
def fasta(p):
    d={};k=None
    for l in p.read_text().splitlines():
        if l.startswith('>'):k=l[1:].split()[0];d[k]=''
        elif l.strip():d[k]+=l.strip()
    return d
core=sorted('LOCAL_'+x['sample'] for x in json.loads((OLD/'genome_sources.json').read_text()) if x['in_core_tree']);assert len(core)==66
genes=['BenA','CaM','RPB2'];data={g:fasta(OLD/'tree_inputs'/(g+'_all_strains.fasta')) for g in genes};manifest=[]
for n in [1,2,3]:
    for combo in itertools.combinations(genes,n):
        name='_'.join(combo);seqs={s:''.join(data[g][s] for g in combo) for s in core};width=len(next(iter(seqs.values())))
        (P/(name+'.fasta')).write_text(''.join('>'+s+'\n'+seqs[s]+'\n' for s in core))
        start=1;parts=[]
        for g in combo:
            end=start+len(data[g][core[0]])-1;parts.append('charset '+g+' = '+str(start)+'-'+str(end)+';');start=end+1
        (P/(name+'.nex')).write_text('#nexus\nbegin sets;\n'+'\n'.join(parts)+'\nend;\n')
        manifest.append({'name':name,'genes':list(combo),'tips':66,'sites':width,'alignment_source':'2026-11-10 validated reference-window alignments; same retained columns; restricted to exact core cohort'})
(P/'design.json').write_text(json.dumps({'core_tips':core,'combinations':manifest,'comparison':'De novo ML for all seven combinations plus pruned-reference-tree sensitivity. Core genome topology is exploratory, not a validated species-tree truth.','rooting':'Unrooted RF; display root SRR10507650 used consistently for the marker combinations, without an ancestral inference.'},indent=2))
print('PREPARED seven core66 alignments')
