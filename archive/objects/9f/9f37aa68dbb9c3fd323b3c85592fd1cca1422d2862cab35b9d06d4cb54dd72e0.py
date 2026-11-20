from pathlib import Path
import sys,os,json,math,csv,gc
OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'comparative_visualization/2026-11-13/vendor'))
from circos_review_common import constrain_process
LIMITS=constrain_process()
os.environ['MPLCONFIGDIR']=str(OUT/'qa/mpl');os.environ['OPENBLAS_NUM_THREADS']='1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from terreus_tree_helpers import parse_newick,support,strong
assert json.loads((OUT/'data/new_tree_validation.json').read_text())['state']=='complete'
root=parse_newick((OUT/'runs/complete230.treefile').read_text());tips=root.leaves();assert len(tips)==230
old=json.loads((ROOT/'comparative_visualization/2026-11-11/analysis/tree_annotations_all_combinations.json').read_text())
ann={r['tip']:dict(r) for r in old if r['tree']=='combined_all_strains' and r['tip'] in {t.name for t in tips}}
palette=json.loads((ROOT/'stakeholder_release/2026-11-17_v3_taxonomy_review/data/palette.json').read_text())['species']
refs={r['reference_id']:r for r in json.loads((ROOT/'taxonomy_reference_build/2026-11-10/analysis/tree_reference_metadata.json').read_text())}
sourcecols={'genome':'#087f8c','reference':'#bd642d'};gold='#e8af27';ink='#213c50'
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','font.size':10,'text.color':ink})
def source(k):return 'genome' if k.startswith('LOCAL_') else 'reference'
def typed(k):return ann[k]['type'] or bool(ann[k].get('genome_type_species'))
def label(k):
 return ('[G] ' if source(k)=='genome' else '[R11] ' if k.startswith('HIST2011_') else '[R] ')+ann[k]['display_label'].replace(' [2011]','').replace('[TG]','[TYPE-G]').replace('[T]','[TYPE]')
# Recalculate type-containing clades on the NEW tree. Old clade IDs are not reused.
leaves={};paths={}
def index(n,path):
 if not n.children:leaves[id(n)]=[n.name];paths[n.name]=path+[n]
 else:
  for c in n.children:index(c,path+[n])
  leaves[id(n)]=sum((leaves[id(c)] for c in n.children),[])
index(root,[])
types={k:r['analysis_species'].replace('Aspergillus_','A. ').replace('_',' ') for k,r in refs.items() if r['type_flag']}
sets={i:{types[k] for k in names if k in types} for i,names in leaves.items()};selected={};chosen={}
for k,path in paths.items():
 pure=[n for n in reversed(path[:-1]) if len(sets[id(n)])==1 and strong(n)]
 mixed=[n for n in reversed(path[:-1]) if 1<len(sets[id(n)])<=8 and strong(n)]
 n=pure[-1] if pure else mixed[0] if mixed else None;chosen[k]=n
 if n is not None:selected[id(n)]=n
ordered=sorted(selected.values(),key=lambda n:(sorted(sets[id(n)]),len(leaves[id(n)])))
codes={id(n):'T'+str(i+1).zfill(2) for i,n in enumerate(ordered)};clades=[]
for n in ordered:
 a,b=support(n);clades.append({'code':codes[id(n)],'type_species':sorted(sets[id(n)]),'SH_aLRT':a,'UFBoot':b,'tips':leaves[id(n)]})
for k,n in chosen.items():
 tt=sets[id(n)] if n is not None else set();a=ann[k]
 a['new_clade_code']=codes[id(n)] if n is not None else '-';a['new_clade_category']=next(iter(tt)) if len(tt)==1 else 'Mixed type species' if tt else 'No supported type clade'
 assert a['available_loci']==3
def cols(k):return [sourcecols[source(k)],'#408f87',palette[ann[k]['species_ring']],palette[ann[k]['new_clade_category']]]
def positions(clado=False):
 ys={id(t):i for i,t in enumerate(tips)};xs={id(root):0}
 def visit(n):
  for c in n.children:xs[id(c)]=xs[id(n)]+(1 if clado else c.length);visit(c)
  if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
 visit(root);return xs,ys
stems=[]
def save(fig,stem):
 if '--legend-only' in sys.argv and stem!='complete230_colour_legend':plt.close(fig);return
 fig.savefig(OUT/'figures'/(stem+'.pdf'),bbox_inches='tight',pad_inches=.25)
 fig.savefig(OUT/'figures'/(stem+'.svg'),bbox_inches='tight',pad_inches=.25)
 fig.savefig(OUT/'qa'/(stem+'.png'),dpi=60,bbox_inches='tight',pad_inches=.25)
 plt.close(fig);gc.collect();stems.append(stem)
fig=plt.figure(figsize=(24,41));ax=fig.add_axes([.035,.03,.34,.91]);lab=fig.add_axes([.395,.03,.38,.91],sharey=ax);track=fig.add_axes([.80,.03,.18,.91],sharey=ax)
xs,ys=positions();mx=max(xs.values())
def branches(n):
 if n.children:
  ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],c='#607384',lw=.7)
  for c in n.children:ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,c='#607384',lw=.7);branches(c)
  if strong(n):ax.scatter(xs[id(n)],ys[id(n)],s=6,c=ink)
 else:
  y=ys[id(n)];ax.plot([xs[id(n)],mx*1.015],[y,y],c='#cbd4dc',lw=.4,ls=':');ax.scatter(mx*1.025,y,marker='s' if source(n.name)=='genome' else 'o',s=15,c=sourcecols[source(n.name)])
branches(root)
for i,t in enumerate(tips):
 k=t.name;lab.text(0,i,label(k),va='center',fontsize=9,color=sourcecols[source(k)],weight='bold' if typed(k) or source(k)=='genome' else 'normal')
 if typed(k):lab.scatter(-.025,i,marker='*',s=55,c=gold,edgecolor='#533f12',linewidth=.4,clip_on=False)
 for j,c in enumerate(cols(k)):track.add_patch(Rectangle((j-.38,i-.45),.76,.90,color=c,lw=0))
 track.text(4.35,i,ann[k]['new_clade_code'],va='center',ha='center',fontsize=8)
for a in [ax,lab,track]:a.set_ylim(len(tips)-.3,-2.2);a.set_yticks([])
ax.set_xlim(0,mx*1.05);ax.spines[['top','right','left']].set_visible(False);ax.set_xlabel('Substitutions per site',fontsize=12)
lab.set_xlim(0,1);lab.axis('off');track.set_xlim(-.55,4.85);track.axis('off')
for x,h in zip([0,1,2,3,4.35],['Sequence\nsource','Three\nloci','Reviewed\nspecies','Reference\ntype clade','Clade\ncode']):track.text(x,1.004,h,transform=track.get_xaxis_transform(),ha='center',va='bottom',fontsize=10)
fig.text(.035,.984,'BenA + CaM + RPB2 | complete-locus cohort',fontsize=27,weight='bold',va='top')
fig.text(.035,.969,'230 tips | 75 genome-derived entries | 155 taxonomic references | 1,000 UFBoot + 1,000 SH-aLRT',fontsize=13,va='top')
save(fig,'complete230_vertical')
fig,ax=plt.subplots(figsize=(38,38),subplot_kw={'projection':'polar'});xs,ys=positions(True);mx=max(xs.values());step=2*math.pi/len(tips)
def circle(n):
 r=xs[id(n)]/mx;a=ys[id(n)]*step
 if n.children:
  aa=[ys[id(c)]*step for c in n.children];ax.plot(np.linspace(min(aa),max(aa),70),[r]*70,c='#607384',lw=.7)
  for c in n.children:ax.plot([ys[id(c)]*step]*2,[r,xs[id(c)]/mx],c='#607384',lw=.7);circle(c)
  if strong(n):ax.scatter(a,r,s=6,c=ink)
 else:
  k=n.name;ax.plot([a]*2,[r,1.025],c='#cbd4dc',lw=.4,ls=':');ax.scatter(a,1.025,marker='s' if source(k)=='genome' else 'o',s=15,c=sourcecols[source(k)])
  for j,c in enumerate(cols(k)):ax.bar(a,.033,width=step*.93,bottom=1.053+j*.045,color=c,lw=0)
  deg=math.degrees(a);left=90<deg<270
  if typed(k):ax.scatter(a,1.245,marker='*',s=70,color=gold,edgecolor='#533f12',linewidth=.4)
  ax.text(a,1.28,label(k),rotation=deg+(180 if left else 0),rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=9,color=sourcecols[source(k)],weight='bold' if typed(k) or source(k)=='genome' else 'normal')
circle(root);ax.set_ylim(0,1.94);ax.axis('off');fig.subplots_adjust(left=.02,right=.98,top=.94,bottom=.01)
fig.text(.035,.983,'BenA + CaM + RPB2 | complete-locus cohort',fontsize=30,weight='bold',va='top');fig.text(.035,.967,'230 tips | circular cladogram | branch lengths not represented',fontsize=14,va='top')
save(fig,'complete230_circular')
fig=plt.figure(figsize=(22,17));la=fig.add_axes([0,0,1,1]);la.set_xlim(0,1);la.set_ylim(0,1);la.axis('off')
def text(x,y,s,size=12,bold=False):return la.text(x,y,s,fontsize=size,weight='bold' if bold else 'normal',va='center')
def swatch(x,y,c,s):la.add_patch(Rectangle((x,y-.005),.009,.010,color=c));text(x+.015,y,s,11)
text(.035,.96,'Complete three-locus tree | shared colour legend',26,True)
text(.035,.915,'Tracks: sequence source  /  three represented loci  /  reviewed species  /  supported reference type clade',14)
swatch(.035,.875,sourcecols['genome'],'[G] Genome-derived');swatch(.30,.875,sourcecols['reference'],'[R] Taxonomic reference; [R11] 2011 reference');swatch(.72,.875,'#408f87','BenA, CaM and RPB2 represented')
la.scatter(.04,.832,marker='*',s=100,c=gold,edgecolor='#533f12');text(.055,.832,'[TYPE] Reference type strain; [TYPE-G] documented genome type',12)
la.scatter(.04,.798,s=20,c=ink);text(.055,.798,'Node dot: SH-aLRT >=80% and UFBoot >=95%; numerical values are retained in the Newick tree',12)
text(.035,.775,'An asterisk after a species name denotes a reported label with unresolved species assignment; gold stars mark types.',11)
text(.035,.753,'Species and type-clade colours',17,True)
names=sorted({ann[k]['species_ring'] for k in ann}|{ann[k]['new_clade_category'] for k in ann})
for i,s in enumerate(names):swatch(.035+(i//12)*.31,.715-(i%12)*.026,palette[s],s)
text(.035,.36,'New-tree type clades (T codes are specific to this tree)',17,True)
for i,c in enumerate(clades):
 x=.035+(i//10)*.31;y=.324-(i%10)*.025
 namestr='; '.join(c['type_species']);namestr=namestr if len(namestr)<43 else 'Mixed types (see clade key)'
 text(x,y,c['code']+'  '+namestr,10)
text(.035,.045,'No missing whole loci or imputation in this cohort. Internal alignment gaps and partial sequence ends remain.',12)
save(fig,'complete230_colour_legend')
(OUT/'data/new_tree_type_clades.json').write_text(json.dumps(clades,indent=2));(OUT/'data/new_tree_annotations.json').write_text(json.dumps(ann,indent=2))
with (OUT/'data/new_tip_colour_audit.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['tip','label','source','available_loci','species','type_clade','type_category','type_marker'])
 for t in tips:w.writerow([t.name,label(t.name),source(t.name),3,ann[t.name]['species_ring'],ann[t.name]['new_clade_code'],ann[t.name]['new_clade_category'],typed(t.name)])
(OUT/'qa/rendered_stems.json').write_text(json.dumps(stems));(OUT/'qa/render_resource_limits.json').write_text(json.dumps(LIMITS));print('Rendered',stems)
