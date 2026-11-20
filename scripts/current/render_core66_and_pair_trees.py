from pathlib import Path
import json,csv,copy,hashlib
import annotate_terreus_trees_v2 as plots
from terreus_tree_helpers import parse_newick,prune,Node,strong
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis';D=P/'core66_combinations';F=P/'figures/core66_trees';F.mkdir(exist_ok=True);PAIR=P/'reference_pair_combinations';keep=set(json.loads((D/'design.json').read_text())['core_tips'])
def newick(n):return ('('+','.join(newick(c) for c in n.children)+')' if n.children else '')+n.name+':'+str(n.length)
def splits(root):
    ss=set()
    def visit(n):
        a=({n.name}&keep) if not n.children else set().union(*(visit(c) for c in n.children));b=keep-a
        if n.children and len(a)>1 and len(b)>1:ss.add(min((tuple(sorted(a)),tuple(sorted(b))),key=lambda x:(len(x),x)))
        return a
    visit(root);return ss
def reroot_outgroup(root,name):
    # Undirected edge representation preserves branch lengths and edge-support labels.
    graph={};objects={}
    def visit(n):
        objects[id(n)]=n;graph.setdefault(id(n),[])
        for c in n.children:
            graph[id(n)].append((id(c),c.length,c.name if c.children else ''));graph.setdefault(id(c),[]).append((id(n),c.length,c.name if c.children else ''));visit(c)
    visit(root);tip=next(id(n) for n in root.leaves() if n.name==name);neighbor,length,edge=graph[tip][0]
    def branch(i,parent,length,label):
        children=[branch(j,i,dist,sup) for j,dist,sup in graph[i] if j!=parent]
        n=Node(objects[i].name if not children else label,length,children)
        if len(children)==1:children[0].length+=n.length;return children[0]
        return n
    return Node('',0,[Node(name,length/2),branch(neighbor,tip,length/2,edge)])
oldann=json.loads((P/'analysis/tree_annotations.json').read_text());allann=list(oldann);allclades=json.loads((P/'analysis/type_clades.json').read_text());summary=[];genome_clades=[];core_annotations=[]
def genome_type_ring(root,ann,tree_name):
    updated=copy.deepcopy(ann);leaves,paths=plots.index_tree(root);anchors={k:plots.annotations[k]['genome_type_species'] for k in keep if plots.annotations[k].get('genome_type_species')};types={i:{anchors[k] for k in ll if k in anchors} for i,ll in leaves.items()};selected={}
    for k,species in anchors.items():
        possible=[n for n in reversed(paths[k][:-1]) if strong(n) and types[id(n)]=={species} and all(updated[x]['type_clade_category']==species or x in anchors and anchors[x]==species for x in leaves[id(n)])]
        if possible:selected[id(possible[-1])]=possible[-1]
        else:genome_clades.append({'tree':tree_name,'status':'type_genome_anchor_only','type_species':species,'tips':[k],'SH_aLRT':None,'UFBoot':None})
    for n in selected.values():
        species=next(iter(types[id(n)]));a,b=plots.support(n);genome_clades.append({'tree':tree_name,'status':'supported_in_tree_type_clade_consistent_with_reference_placement','type_species':species,'tips':leaves[id(n)],'SH_aLRT':a,'UFBoot':b})
        for k in leaves[id(n)]:updated[k]['genome_type_clade_category']=species
    for k in keep:core_annotations.append({**updated[k],'tree':tree_name})
    return updated
genome=parse_newick(Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\0_4_mafft_iqtree\Sep28\concatenated.fa.treefile').read_text())
for t in genome.leaves():t.name='LOCAL_'+t.name.removesuffix('_transcripts').removesuffix('.transcripts')
target=splits(genome)
for d in json.loads((D/'design.json').read_text())['combinations']:
    name=d['name'];source_name='combined_all_strains' if len(d['genes'])==3 else name
    refpath=(PAIR/(name+'.treefile')) if len(d['genes'])==2 else OLD/'trees'/(source_name+'.treefile')
    iq=refpath.with_suffix('.iqtree')
    if not refpath.exists() or not iq.exists() or 'Date and time:' not in iq.read_text():print('WAITING_REFERENCE',name,flush=True);continue
    reference=parse_newick(refpath.read_text());ann,clades=plots.type_clades(reference,source_name)
    if len(d['genes'])==2:
        allann+=list(ann.values());allclades+=clades;plots.F=P/'figures/trees'
        stamp=plots.F/(name+'.source_sha256.txt');digest='genome_types_v2:'+hashlib.sha256(refpath.read_bytes()).hexdigest()
        if not stamp.exists() or stamp.read_text()!=digest:
            for focused in [False,True]:
                allowed={k for k,a in plots.annotations.items() if k.startswith('LOCAL_') or a['type'] or (k in plots.meta and (plots.meta[k].get('accepted_type_name') or plots.meta[k].get('historical_type_name')))}
                r=prune(reference,allowed) if focused else reference;stem=name+('_type_focused' if focused else '_all_references');title=name.replace('_',' + ')+(' | type-focused display' if focused else ' | all reference records');plots.linear(r,title,stem,ann,clades);plots.circular(r,title,stem,ann,clades)
            stamp.write_text(digest)
    reduced=prune(reference,keep);de_novo=parse_newick((D/(name+'.treefile')).read_text());assert {t.name for t in reduced.leaves()}==keep
    (D/(name+'.reference_pruned66.treefile')).write_text(newick(reduced)+';\n');plots.F=F
    for source,root in [('de_novo',de_novo),('reference_pruned',reduced)]:
        before=splits(root);root=reroot_outgroup(root,'LOCAL_SRR10507650');assert splits(root)==before;stem=name+'_'+source
        localann=genome_type_ring(root,ann,stem);stamp=F/(stem+'.genome_types_v2.done')
        if not stamp.exists():
            title=name.replace('_',' + ')+(' | inferred on 66 genomes' if source=='de_novo' else ' | reference tree reduced to 66 genomes');plots.linear(root,title,stem,localann,[],True);plots.circular(root,title,stem,localann,[],True);stamp.write_text('completed')
    s1,s2=splits(de_novo),splits(reduced);summary.append({'combination':name,'shared_tips':66,'de_novo_RF_to_core':len(s1^target)/126,'reference_pruned_RF_to_core':len(s2^target)/126,'de_novo_vs_reference_pruned_RF':len(s1^s2)/126,'reference_tree':str(refpath),'reference_tree_sha256':hashlib.sha256(refpath.read_bytes()).hexdigest()});print('RENDERED CORE66',name,flush=True)
ann={x['tip']:x for x in oldann if x['tree']=='combined_all_strains'};plots.F=F
root=reroot_outgroup(genome,'LOCAL_SRR10507650');localann=genome_type_ring(root,ann,'core_genome')
if not (F/'core_genome.genome_types_v2.done').exists():
    plots.linear(root,'Core genome | original 4,661-locus topology','core_genome',localann,[],True);plots.circular(root,'Core genome | original 4,661-locus topology','core_genome',localann,[],True);(F/'core_genome.genome_types_v2.done').write_text('completed')
(P/'analysis/tree_annotations_all_combinations.json').write_text(json.dumps(allann,indent=2));(P/'analysis/type_clades_all_combinations.json').write_text(json.dumps(allclades,indent=2));(D/'reference_sampling_sensitivity.json').write_text(json.dumps(summary,indent=2))
(P/'analysis/core66_genome_type_clades.json').write_text(json.dumps(genome_clades,indent=2));(P/'analysis/core66_tree_annotations.json').write_text(json.dumps(core_annotations,indent=2))
for fn,rows in [('tree_annotations_all_combinations.tsv',allann),('type_clades_all_combinations.tsv',[{**x,'descendant_tips':';'.join(x['descendant_tips'])} for x in allclades])]:
    with (P/'analysis'/fn).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('READY',len(summary),'/7 combinations')
