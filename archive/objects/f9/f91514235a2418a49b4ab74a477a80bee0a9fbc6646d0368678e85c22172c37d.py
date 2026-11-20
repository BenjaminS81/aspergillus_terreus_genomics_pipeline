"""Type-exclusive monophyly annotation, with topology and support kept separate."""
from pathlib import Path
import os,sys,json,hashlib,math,collections,copy,gc,textwrap
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
from terreus_tree_helpers import Node,parse_newick,support,strong,prune
OUT=ROOT/'stakeholder_release/2026-11-20_v1_type_monophyly'
for d in ['data','figures','qa','scripts']:(OUT/d).mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(OUT/'qa/mpl')
PREV=ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene'
CORE=ROOT/'stakeholder_release/2026-11-17_v1/data/core_genome_original.treefile'
GENE=PREV/'runs/complete230.treefile'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(name,x):(OUT/'data'/name).write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def short(s):return s.replace('Aspergillus_','A. ').replace('_',' ')
ann=read(PREV/'data/new_tree_annotations.json');refs={r['reference_id']:r for r in read(ROOT/'taxonomy_reference_build/2026-11-10/analysis/tree_reference_metadata.json')}
gene=parse_newick(GENE.read_text());core0=parse_newick(CORE.read_text())
core_names={t.name.removesuffix('_transcripts') for t in core0.leaves()}
for t in core0.leaves():t.name='LOCAL_'+t.name.removesuffix('_transcripts')
def idx(root):
 nodes=[];leaves={};parents={}
 def visit(n,p=None):
  nodes.append(n);parents[id(n)]=p
  for c in n.children:visit(c,n)
  leaves[id(n)]={x.name for x in n.leaves()}
 visit(root);return nodes,leaves,parents
gn,gl,gp=idx(gene);cn,cl,cp=idx(core0)
types={k:short(r['analysis_species']) for k,r in refs.items() if r['type_flag'] and k in gl[id(gene)]}
for k,a in ann.items():
 if a.get('genome_type_species'):types[k]=a['genome_type_species']
assert len(types)==32
groups=[];assign={}
for species in sorted(set(types.values())):
 anchors={k for k,v in types.items() if v==species}
 candidates=[n for n in gn if anchors<=gl[id(n)] and {types[k] for k in gl[id(n)] if k in types}=={species}]
 if not candidates:
  groups.append({'species':species,'status':'conflicting_type_anchors','anchors':sorted(anchors),'tips':[]});continue
 n=max(candidates,key=lambda n:len(gl[id(n)]));ls=gl[id(n)];a,b=support(n)
 status='supported' if n.children and strong(n) else 'topology_only' if n.children else 'type_tip_only'
 conflicts=[{'tip':k,'reported_reference_species':short(refs[k]['analysis_species'])} for k in sorted(ls) if k in refs and short(refs[k]['analysis_species'])!=species]
 group={'species':species,'status':status,'anchors':sorted(anchors),'anchor_cultures':[ann[k].get('genome_type_culture') or ann[k]['strain'] for k in sorted(anchors)],'tips':sorted(ls),'genome_tips':sorted(k for k in ls if k.startswith('LOCAL_')),'SH_aLRT':a,'UFBoot':b,'conflicting_non_type_reference_labels':conflicts,'node':n}
 groups.append(group)
 for k in ls:
  assert k not in assign,(k,species)
  assign[k]=group
# All matched core tips must be represented by actual loci in the current reference tree.
assert {'LOCAL_'+s for s in core_names}<=gl[id(gene)]
allcore=cl[id(core0)]
def split_key(s,alltips):
 other=alltips-s
 return tuple(sorted(s)) if (len(s),tuple(sorted(s)))<=(len(other),tuple(sorted(other))) else tuple(sorted(other))
def splits(root):
 ns,ls,_=idx(root);alltips=ls[id(root)]
 return {split_key(ls[id(n)],alltips) for n in ns if n is not root and 1<len(ls[id(n)])<len(alltips)-1}
for g in groups:
 members=set(g['tips'])&allcore
 g['core_tips']=sorted(members)
 g['core_unrooted_monophyly']=None if not members else True if len(members)==1 else any(members==s or members==allcore-s for s in cl.values())
 if members:
  containing=min((n for n in cn if members<=cl[id(n)]),key=lambda n:len(cl[id(n)]))
  g['core_original_root_mrca_extra_tips']=sorted(cl[id(containing)]-members)
# Orient the core tree for display using the deepest matched split in the reference topology.
matched=prune(gene,allcore);base_sets=[{t.name for t in c.leaves()} for c in matched.children]
display_outgroup=min(base_sets,key=len)
edge=next((n for n in cn if n is not core0 and (cl[id(n)]==display_outgroup or allcore-cl[id(n)]==display_outgroup)),None)
assert edge is not None,'Induced reference display split absent from core; review rooting rather than forcing it'
adj={id(n):[] for n in cn};lookup={id(n):n for n in cn}
for n in cn:
 for c in n.children:
  adj[id(n)].append((id(c),c.length,c.name if c.children else ''));adj[id(c)].append((id(n),c.length,c.name if c.children else ''))
u=id(edge);v=id(cp[u]);edge_length=edge.length
def orient(cur,prev,length,label):
 old=lookup[cur];children=[orient(other,cur,l,lab) for other,l,lab in adj[cur] if other!=prev]
 return Node(old.name if not old.children else label,length,children)
core=Node('',0,[orient(u,v,edge_length/2,edge.name),orient(v,u,edge_length/2,edge.name)])
# Collapse any degree-two original root without affecting distances.
def simplify(n):
 n.children=[simplify(c) for c in n.children]
 if len(n.children)==1:
  c=n.children[0];c.length+=n.length;return c
 return n
core=simplify(core)
assert splits(core)==splits(core0)
def pair_distances(root):
 paths={}
 def walk(n,path,depth):
  path=path+[(id(n),depth)]
  if not n.children:paths[n.name]=path
  for c in n.children:walk(c,path,depth+c.length)
 walk(root,[],0);result={};names=sorted(paths)
 for i,a in enumerate(names):
  da=dict(paths[a])
  for b in names[i+1:]:
   common=next(depth for node,depth in reversed(paths[b]) if node in da)
   result[(a,b)]=paths[a][-1][1]+paths[b][-1][1]-2*common
 return result
original_distances=pair_distances(core0);display_distances=pair_distances(core)
max_distance_change=max(abs(v-display_distances[k]) for k,v in original_distances.items())
assert max_distance_change<1e-9,max_distance_change
def order(n):
 for c in n.children:order(c)
 n.children.sort(key=lambda c:(-len(c.leaves()),min(t.name for t in c.leaves())))
order(core);order(gene)
def newick(n,rename=False):
 label=n.name
 if not n.children and rename:
  g=assign.get(label);species=g['species'] if g else 'Unresolved'
  label=species.replace('A. ','Aspergillus_').replace(' ','_')+'__'+label.removeprefix('LOCAL_')
 return ('('+','.join(newick(c,rename) for c in n.children)+')' if n.children else '')+label+':'+format(n.length,'.12g')
(OUT/'data/core66_display_rooted.treefile').write_text(newick(core)+';')
(OUT/'data/core66_type_annotated.treefile').write_text(newick(core,True)+';')
(OUT/'data/three_gene230_type_annotated.treefile').write_text(newick(gene,True)+';')
audit=[]
for k in sorted(gl[id(gene)]):
 g=assign.get(k);a=ann[k];species=g['species'] if g else 'Unresolved'
 audit.append({'tip':k,'sample':k.removeprefix('LOCAL_'),'source':'Genome' if k.startswith('LOCAL_') else 'Reference','reported_species':a['taxon'],'assigned_species':species,'confidence':g['status'] if g else 'unresolved','type_anchor_tips':g['anchors'] if g else [],'type_cultures':g['anchor_cultures'] if g else [],'SH_aLRT':g.get('SH_aLRT') if g else None,'UFBoot':g.get('UFBoot') if g else None,'is_type':k in types,'in_core66':k in allcore,'core_monophyly':g.get('core_unrooted_monophyly') if g else None,'reported_name_changed':species!='Unresolved' and species!=a['taxon']})
counts=collections.Counter((assign[k]['species'],assign[k]['status']) if k in assign else ('Unresolved','unresolved') for k in allcore)
summary={'core_tips':66,'three_gene_tips':230,'genome_tips':75,'reference_tips':155,'type_anchors':len(types),'core_assignments':[{'species':k[0],'confidence':k[1],'count':v} for k,v in sorted(counts.items())],'core_reported_names_changed':sum(r['in_core66'] and r['reported_name_changed'] for r in audit),'display_root_split':sorted(display_outgroup),'unrooted_core_topology_preserved':True,'rule':'For each species, select the largest descendant clade containing every included verified type anchor for that species and no type anchor of another species. Include every descendant query/reference tip. Non-type reported names do not define boundaries; conflicts are explicitly retained. Topological membership and SH-aLRT/UFBoot support are separate fields. A type tip alone does not classify neighbouring tips.','support_rule':'supported = SH-aLRT >=80 and UFBoot >=95 at the enclosing clade; topology_only otherwise. These are analysis thresholds, not formal nomenclatural proof.','core_transfer_rule':'Exact existing genome sample mapping only. Transfer the three-gene type-derived name to the corresponding core tip; do not propagate across a mixed core clade. Test each transferred group against unrooted core splits.','rooting_note':'Core rerooted only for display on the deepest split induced by the three-gene tree among matched tips. This is not an independent ancestral-root validation. Original unrooted topology and branch lengths retained.','new_downloads_included':False,'resource_limits':LIMITS,'inputs':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [CORE,GENE,PREV/'data/new_tree_annotations.json',ROOT/'taxonomy_reference_build/2026-11-10/analysis/tree_reference_metadata.json']]}
save('assignment_audit.json',audit);save('type_monophyly_clades.json',[{k:v for k,v in g.items() if k!='node'} for g in groups]);save('analysis_summary.json',summary)
save('rooting_validation.json',{'all_pairwise_distances_checked':len(original_distances),'maximum_absolute_distance_change':max_distance_change,'unrooted_splits_identical':True})
print(json.dumps(summary['core_assignments'],indent=2));print('Core topology conflicts',[(g['species'],len(g['core_tips'])) for g in groups if g.get('core_unrooted_monophyly') is False]);print('Changed',summary['core_reported_names_changed'])
if '--audit-only' in sys.argv:sys.exit()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
INK='#203849';GRAY='#70818f';GOLD='#d9a72f'
palette=read(ROOT/'stakeholder_release/2026-11-17_v3_taxonomy_review/data/palette.json')['species']
# Keep named-species hues consistent with the preceding release; use dark text for readability.
sourcecols={'Genome':'#087f8c','Reference':'#bd642d'}
confidencecols={'supported':'#244f67','topology_only':'#d39a35','type_tip_only':'#9b80ad','unresolved':'#bfc5ca'}
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','font.size':10,'text.color':INK,'axes.labelcolor':INK})
auditmap={r['tip']:r for r in audit};outputs=[]
def species(k):return assign[k]['species'] if k in assign else 'Unresolved'
def color(k):return palette[species(k)]
def label(k,source=False):
 r=auditmap[k];a=ann[k];name=r['assigned_species']
 if name=='Unresolved':name='Aspergillus sp.'
 suffix=' [T-G]' if a.get('genome_type_species') else ' [T]' if k in types else ''
 return ('[G] ' if source and k.startswith('LOCAL_') else '[R] ' if source else '')+name+' | '+a['strain']+suffix
def positions(root,clado=False):
 tips=root.leaves();ys={id(t):i for i,t in enumerate(tips)};xs={id(root):0}
 def walk(n):
  for c in n.children:xs[id(c)]=xs[id(n)]+(1 if clado else c.length);walk(c)
  if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
 walk(root);return tips,xs,ys
def savefig(fig,stem):
 fig.savefig(OUT/'figures'/(stem+'.pdf'),facecolor='white')
 fig.savefig(OUT/'figures'/(stem+'.svg'),facecolor='white')
 plt.close(fig);gc.collect();outputs.append(stem);print('PDF/SVG',stem,flush=True)
def pure_blocks(root):
 ns,ls,pa=idx(root);found=[]
 for n in ns:
  names={species(k) for k in ls[id(n)]}
  if len(names)==1 and 'Unresolved' not in names:
   par=pa[id(n)];parentnames={species(k) for k in ls[id(par)]} if par else set()
   if par is None or len(parentnames)>1:found.append((n,next(iter(names))))
 return found
def linear(root,which):
 tips,xs,ys=positions(root);N=len(tips);isgene=which=='three_gene';mx=max(xs.values())
 fig=plt.figure(figsize=(19,34) if isgene else (18,15))
 ax=fig.add_axes([.035,.055,.345,.88]);lab=fig.add_axes([.395,.055,.355,.88],sharey=ax);tr=fig.add_axes([.755,.055,.055,.88],sharey=ax);br=fig.add_axes([.83,.055,.16,.88],sharey=ax)
 def draw(n):
  if n.children:
   ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],c=GRAY,lw=.65)
   for c in n.children:
    names={species(t.name) for t in c.leaves()};co=palette[next(iter(names))] if len(names)==1 else GRAY
    ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,c=co,lw=.85);draw(c)
   if isgene and strong(n):ax.scatter(xs[id(n)],ys[id(n)],s=5,c=INK,zorder=4)
  else:ax.plot([xs[id(n)],mx*1.025],[ys[id(n)]]*2,c='#dce1e5',lw=.4,ls=':')
 draw(root)
 for i,t in enumerate(tips):
  k=t.name;r=auditmap[k];g=assign.get(k);fs=8.2 if isgene else 10
  lab.text(.025,i,label(k,isgene),fontsize=fs,va='center',weight='bold' if k in types else 'normal',color=INK)
  lab.scatter(.006,i,marker='s' if k.startswith('LOCAL_') else 'o',s=10 if isgene else 18,c=sourcecols[r['source']],clip_on=False)
  for j,co in enumerate([color(k),confidencecols[r['confidence']]]):tr.add_patch(Rectangle((j+.08,i-.44),.72,.88,color=co,lw=0))
  if not isgene:tr.add_patch(Rectangle((2.08,i-.44),.72,.88,color='#b76b40' if r['core_monophyly'] is False else '#4b777c',lw=0))
  if k in types:lab.scatter(-.013,i,marker='*',s=38 if isgene else 65,c=GOLD,edgecolor='#674c16',lw=.35,clip_on=False)
 blocks=[(g['node'],g['species']) for g in groups if g.get('node') and len(g['tips'])>1] if isgene else pure_blocks(root)
 ns,ls,_=idx(root);tipy={t.name:i for i,t in enumerate(tips)}
 for n,sp in blocks:
  vals=[tipy[k] for k in ls[id(n)]];lo=min(vals);hi=max(vals);mid=(lo+hi)/2
  br.plot([.05,.0,.0,.05],[lo-.35,lo-.35,hi+.35,hi+.35],c=palette[sp],lw=2)
  g=next(g for g in groups if g['species']==sp);weak=g['status']=='topology_only'
  text=sp+(' *' if weak else '')
  br.text(.075,mid,text,fontsize=12 if isgene else 12.4,weight='bold',fontstyle='italic',va='center',color=INK)
 for a in [ax,lab,tr,br]:a.set_ylim(N-.35,-1.7);a.set_yticks([])
 ax.set_xlim(0,mx*1.04);ax.spines[['top','right','left']].set_visible(False);ax.set_xlabel('Substitutions per site',fontsize=11)
 lab.set_xlim(0,1);lab.axis('off');tr.set_xlim(0,2 if isgene else 3);tr.axis('off');br.set_xlim(0,1);br.axis('off')
 for j,title in enumerate(['Species','Support']+([] if isgene else ['Core'])):tr.text(j+.44,1.006,title,transform=tr.get_xaxis_transform(),rotation=90,va='bottom',ha='center',fontsize=8)
 title='BenA + CaM + RPB2' if isgene else 'Core-genome phylogeny'
 fig.text(.035,.98,title+' | type-anchored species',fontsize=23,weight='bold',va='top')
 fig.text(.035,.962,('230 tips | 75 genome-derived entries | 155 taxonomic references' if isgene else '66 genomes | species names transferred from the complete three-gene tree'),fontsize=11.5,va='top')
 savefig(fig,'three_gene230_type_monophyly_linear' if isgene else 'core66_type_monophyly_linear')
def circular(root,which):
 tips,xs,ys=positions(root,True);N=len(tips);isgene=which=='three_gene';mx=max(xs.values());step=2*math.pi/N
 fig,ax=plt.subplots(figsize=(29,29) if isgene else (21,21),subplot_kw={'projection':'polar'})
 def draw(n):
  rad=xs[id(n)]/mx;angle=ys[id(n)]*step
  if n.children:
   angles=[ys[id(c)]*step for c in n.children];ax.plot(np.linspace(min(angles),max(angles),50),[rad]*50,c=GRAY,lw=.65)
   for c in n.children:ax.plot([ys[id(c)]*step]*2,[rad,xs[id(c)]/mx],c=GRAY,lw=.65);draw(c)
  else:
   k=n.name;r=auditmap[k];ax.plot([angle]*2,[rad,1.015],c='#cbd4dc',lw=.4,ls=':')
   colors=[sourcecols[r['source']],color(k),confidencecols[r['confidence']]]+([] if isgene else ['#b76b40' if r['core_monophyly'] is False else '#4b777c'])
   for j,co in enumerate(colors):ax.bar(angle,.028,width=step*.94,bottom=1.035+j*.041,color=co,lw=0)
   labelrad=1.21 if isgene else 1.25
   if k in types:ax.scatter(angle,labelrad-.022,marker='*',s=42,c=GOLD,edgecolor='#674c16',lw=.3)
   deg=math.degrees(angle);left=90<deg<270
   ax.text(angle,labelrad,label(k,isgene),rotation=deg+180 if left else deg,rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=7.8 if isgene else 9.2,weight='bold' if k in types else 'normal')
 draw(root);ax.set_ylim(0,1.9 if isgene else 2.0);ax.axis('off');fig.subplots_adjust(left=.015,right=.985,bottom=.02,top=.93)
 fig.text(.035,.98,('BenA + CaM + RPB2' if isgene else 'Core-genome phylogeny')+' | type-anchored species',fontsize=25,weight='bold',va='top')
 fig.text(.035,.959,f'{N} tips | circular cladogram | branch lengths not represented',fontsize=12,va='top')
 savefig(fig,('three_gene230' if isgene else 'core66')+'_type_monophyly_circular')
linear(core,'core');circular(core,'core');linear(gene,'three_gene');circular(gene,'three_gene')
# Matched-tip comparison: no references or missing genomes introduced into either topology.
reduced=prune(gene,allcore);order(reduced)
fig=plt.figure(figsize=(19,17));left=fig.add_axes([.035,.05,.25,.86]);middle=fig.add_axes([.30,.05,.40,.86]);right=fig.add_axes([.715,.05,.25,.86])
orders=[]
for root,ax,reverse in [(core,left,False),(reduced,right,True)]:
 tips,xs,ys=positions(root,True);mx=max(xs.values());orders.append({t.name:i for i,t in enumerate(tips)})
 def draw(n):
  if n.children:
   ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],c=GRAY,lw=.65)
   for c in n.children:ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,c=GRAY,lw=.65);draw(c)
  else:ax.plot([xs[id(n)],mx+.5],[ys[id(n)]]*2,c=color(n.name),lw=.8)
 draw(root);ax.set_ylim(65.7,-1);ax.set_xlim(mx+.6,0) if reverse else ax.set_xlim(0,mx+.6);ax.axis('off')
for k,y in orders[0].items():
 yy=orders[1][k];middle.plot([.18,.82],[y,yy],c=color(k),lw=.65,alpha=.55)
 middle.text(.16,y,k.removeprefix('LOCAL_'),ha='right',va='center',fontsize=8);middle.text(.84,yy,k.removeprefix('LOCAL_'),ha='left',va='center',fontsize=8)
middle.set_xlim(0,1);middle.set_ylim(65.7,-1);middle.axis('off')
fig.text(.035,.975,'Core genome and three-gene topology | identical 66-genome set',fontsize=22,weight='bold',va='top');fig.text(.035,.945,'Core genome',fontsize=13,weight='bold');fig.text(.76,.945,'BenA + CaM + RPB2',fontsize=13,weight='bold')
savefig(fig,'matched66_type_monophyly_comparison')
# A single shared legend, including the meaning of weak labels and core discordance.
fig=plt.figure(figsize=(17,12));ax=fig.add_axes([0,0,1,1]);ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
ax.text(.04,.95,'Shared species colours and evidence key',fontsize=23,weight='bold',va='top')
names=sorted({species(k) for k in gl[id(gene)]})
for i,sp in enumerate(names):
 x=.04+(i//15)*.35;y=.87-(i%15)*.037
 ax.add_patch(Rectangle((x,y-.008),.014,.016,color=palette[sp]));ax.text(x+.024,y,sp,va='center',fontsize=11)
entries=[('Genome-derived [G]','#087f8c'),('Taxonomic reference [R]','#bd642d'),('Supported type clade','#244f67'),('Topology-only type clade *','#d39a35'),('Type tip alone','#9b80ad'),('Unassigned','#bfc5ca'),('Core group monophyletic / singleton','#4b777c'),('Core group not monophyletic','#b76b40')]
for i,(txt,co) in enumerate(entries):
 x=.74;y=.87-i*.045;ax.add_patch(Rectangle((x,y-.008),.014,.016,color=co));ax.text(x+.024,y,txt,va='center',fontsize=10)
ax.scatter(.045,.265,marker='*',s=110,c=GOLD,edgecolor='#674c16');ax.text(.065,.265,'[T] taxonomic reference type; [T-G] documented genome type',fontsize=11,va='center')
lines=['Strict type-exclusive rule: include all descendants of the largest clade containing the species type anchor(s) and no other species type.',
 'Support is separate from topology: supported = SH-aLRT >=80 and UFBoot >=95; * marks a type-exclusive clade below either threshold.',
 'Core names are transferred by matched genome ID. Brown Core cells identify groups that do not form a single unrooted core split.',
 'Linear tracks: species / clade support / Core (core tree only). Circular tracks, inner to outer: source / species / support / Core.',
 'Right-hand names are centred on clades. Core names may repeat across separate pure blocks; mixed clades are not renamed wholesale.',
 'Node dots on the three-gene linear tree mark SH-aLRT >=80 and UFBoot >=95. Circular trees are cladograms; linear trees retain branch lengths.',
 'Display rooting follows the deepest reference-tree split among the matched genomes; it does not validate an ancestral root.',
 'Names are phylogenetic working assignments. Reported labels, type evidence and non-type reference conflicts remain in the audit.',
 'Assignments are conditional on sampled types: A. jilinensis and A. ralucae types lack complete three-locus data and are not in this tree.']
for i,line in enumerate(lines):ax.text(.04,.225-i*.024,line,fontsize=10,va='top')
savefig(fig,'shared_species_and_evidence_legend')
save('figure_stems.json',outputs)
print('FIGURES COMPLETE',flush=True)
