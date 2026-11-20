from pathlib import Path
import json,csv,collections
from terreus_tree_helpers import *
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
meta={d['reference_id']:d for d in json.loads((P/'tree_reference_metadata.json').read_text())};types={k:d['analysis_species'] for k,d in meta.items() if d['type_flag']};calls=json.loads((P/'genomic_markers/locus_calls.json').read_text());all_samples=sorted({d['sample'] for d in calls});results=[]
for tree in ['combined_all_strains','BenA','CaM','RPB2','type_focused_sensitivity']:
    assert (P/'trees'/(tree+'.iqtree')).exists(),tree+' not completed'
    root=parse_newick((P/'trees'/(tree+'.treefile')).read_text());paths={};type_sets={};leaves={}
    def visit(n,path):
        if not n.children:paths[n.name]=path+[n];leaves[id(n)]=[n.name]
        else:
            for c in n.children:visit(c,path+[n])
            leaves[id(n)]=[x for c in n.children for x in leaves[id(c)]]
        type_sets[id(n)]={types[x] for x in leaves[id(n)] if x in types}
    visit(root,[])
    for ident,path in paths.items():
        if not ident.startswith('LOCAL_'):continue
        pure=[n for n in reversed(path[:-1]) if len(type_sets[id(n)])==1 and strong(n)]
        broad=[n for n in reversed(path[:-1]) if 0<len(type_sets[id(n)])<=8 and strong(n)]
        node=pure[-1] if pure else broad[0] if broad else None
        results.append({'sample':ident[6:],'tree':tree,'status':'supported_single_type_species_clade' if pure else 'supported_multiple_type_species_clade' if node else 'no_resolved_supported_type_clade','type_species':','.join(sorted(type_sets[id(node)])) if node else '', 'SH_aLRT':support(node)[0] if node else '', 'UFBoot':support(node)[1] if node else '', 'descendant_tips':len(leaves[id(node)]) if node else ''})
(P/'supported_type_clades.json').write_text(json.dumps(results,indent=2))
rows=[]
for sample in all_samples:
    row={'sample':sample};unique=[];available=0
    for tree in ['combined_all_strains','BenA','CaM','RPB2','type_focused_sensitivity']:
        d=next((d for d in results if d['sample']==sample and d['tree']==tree),{})
        for field in ['status','type_species','SH_aLRT','UFBoot']:row[tree+'_'+field]=d.get(field,'not_in_tree' if field=='status' else '')
        if tree in ['BenA','CaM','RPB2'] and d:
            available+=1
            if d['status']=='supported_single_type_species_clade':unique.append(d['type_species'])
    row['available_loci']=available
    row['interpretation']='no_validated_locus' if available==0 else 'single_locus_only' if available==1 else 'supported_locus_conflict_review' if len(set(unique))>1 else 'three_locus_supported_concordance' if len(unique)==3 else 'incomplete_species_resolution_review'
    row['supported_single_species_loci']=len(unique);row['concordant_supported_species']=unique[0] if len(set(unique))==1 else ''
    full=row['combined_all_strains_type_species'];reduced=row['type_focused_sensitivity_type_species']
    row['full_vs_sensitivity']='same_supported_type_clade' if full and full==reduced else 'supported_single_species_conflict' if row['combined_all_strains_status']==row['type_focused_sensitivity_status']=='supported_single_type_species_clade' and full!=reduced else 'support_or_resolution_differs'
    rows.append(row)
(P/'sample_taxonomic_review.json').write_text(json.dumps(rows,indent=2))
with (P/'sample_taxonomic_review.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('TAXONOMIC STATUS',collections.Counter(d['interpretation'] for d in rows))
for d in rows:
    if d['sample'] in ['IFO6365','TN-484','AT2022','SRR25636141','SRR25636143','NIH2624','HNGD-TM15']:print(d)
