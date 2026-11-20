"""Correct the core display without altering taxon assignments or tree topology."""
from pathlib import Path
import sys,json,hashlib,collections,gc,shutil,os
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
limits=constrain_process()
from terreus_tree_helpers import parse_newick
OLD=ROOT/'stakeholder_release/2026-11-20_v1_type_monophyly'
OUT=ROOT/'stakeholder_release/2026-11-20_v2_core_label_fix'
for d in ['figures','data','qa','scripts']:(OUT/d).mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(OLD/'qa/mpl')
def load(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
treepath=OLD/'data/core66_display_rooted.treefile';r=parse_newick(treepath.read_text())
audit=load(OLD/'data/assignment_audit.json');ann={x['tip']:x for x in audit}
groups=load(OLD/'data/type_monophyly_clades.json');palette=load(ROOT/'stakeholder_release/2026-11-17_v3_taxonomy_review/data/palette.json')['species']
tips=r.leaves();names={t.name for t in tips};assert len(tips)==len(names)==66
nodes=[];leaves={};xs={id(r):0};ys={id(t):i for i,t in enumerate(tips)}
def index(n):
 nodes.append(n);leaves[id(n)]={t.name for t in n.leaves()}
 for c in n.children:xs[id(c)]=xs[id(n)]+c.length;index(c)
 if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
index(r)
side=[];status=[];bad=set()
for sp in sorted({ann[k]['assigned_species'] for k in names}):
 target={k for k in names if ann[k]['assigned_species']==sp}
 mono=any(target==s or target==names-s for s in leaves.values())
 if len(target)>1 and mono:
  node=next((n for n in nodes if leaves[id(n)]==target),None)
  assert node is not None,'Display root bisects a pure group; review rather than label a non-clade'
  side.append((node,sp,palette[sp]));state='single monophyletic core group'
 elif len(target)==1:state='singleton; no multi-tip monophyly claim'
 else:state='not monophyletic in current core tree';bad.update(target)
 status.append({'species':sp,'count':len(target),'core_status':state})
mixed=min((n for n in nodes if bad<=leaves[id(n)]),key=lambda n:len(leaves[id(n)]))
assert leaves[id(mixed)]==bad and len(bad)==40
assert {ann[k]['assigned_species'] for k in bad}=={'A. terreus','A. montoensis'}
side.append((mixed,'Mixed three-gene\nassignments\nterreus / montoensis','#687985'))
gene=parse_newick((ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene/runs/complete230.treefile').read_text());gene_nodes=[]
def scan(n):
 gene_nodes.append(n)
 for c in n.children:scan(c)
scan(gene)
for g in groups:
 if g['tips']:assert any({t.name for t in n.leaves()}==set(g['tips']) for n in gene_nodes),g['species']
record={'correction':'Previously every maximal same-label core subgroup received a large species label. The correction labels only globally monophyletic multi-tip core groups; the 40-tip intermingled terreus/montoensis region gets one neutral mixed-assignment bracket.','previous_repeated_side_labels':{'A. terreus':11,'A. montoensis':8},'tip_labels_unchanged':True,'unique_core_tips':66,'three_gene_clade_memberships_rechecked':True,'tree_source':str(treepath),'tree_sha256':hashlib.sha256(treepath.read_bytes()).hexdigest(),'source_assignment_sha256':hashlib.sha256((OLD/'data/assignment_audit.json').read_bytes()).hexdigest(),'core_groups':status,'large_side_labels':[label for _,label,_ in side],'species_name_asterisks_removed_from_core_display':True,'support_encoding':'Dedicated three-gene support column. No large core species-name asterisks. Gold leaf stars mark documented type genomes only.','limits':limits,'interpretation':'This is discordance between the current stored topology and transferred assignments, not proof of biological polyphyly. Sequence provenance, locus orthology and core-tree validation remain separate gates.'}
save(OUT/'data/correction_audit.json',record);save(OUT/'data/core_tip_assignments.json',[a for a in audit if a['in_core66']])
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
INK='#203849';GRAY='#70818f';GOLD='#d9a72f';SUP={'supported':'#244f67','topology_only':'#d39a35'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','text.color':INK})
def output(fig,stem):
 fig.savefig(OUT/'figures'/(stem+'.pdf'),facecolor='white');fig.savefig(OUT/'figures'/(stem+'.svg'),facecolor='white');plt.close(fig);gc.collect()
fig=plt.figure(figsize=(18,15));ax=fig.add_axes([.035,.055,.345,.88]);lab=fig.add_axes([.395,.055,.355,.88],sharey=ax);tr=fig.add_axes([.755,.055,.055,.88],sharey=ax);br=fig.add_axes([.83,.055,.16,.88],sharey=ax);mx=max(xs.values())
def branches(n):
 if n.children:
  ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],c=GRAY,lw=.65)
  for c in n.children:
   ss={ann[k]['assigned_species'] for k in leaves[id(c)]};co=palette[next(iter(ss))] if len(ss)==1 else GRAY
   ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,c=co,lw=.85);branches(c)
 else:ax.plot([xs[id(n)],mx*1.025],[ys[id(n)]]*2,c='#dce1e5',lw=.4,ls=':')
branches(r)
for i,t in enumerate(tips):
 a=ann[t.name];label=a['assigned_species']+' | '+a['sample']+(' [T-G]' if a['is_type'] else '')
 lab.text(.025,i,label,fontsize=10,va='center',weight='bold' if a['is_type'] else 'normal',color=INK);lab.scatter(.006,i,marker='s',s=18,c='#087f8c')
 if a['is_type']:lab.scatter(-.013,i,marker='*',s=65,c=GOLD,edgecolor='#674c16',lw=.35,clip_on=False)
 colors=[palette[a['assigned_species']],SUP[a['confidence']],'#b76b40' if a['core_monophyly'] is False else '#4b777c']
 for j,co in enumerate(colors):tr.add_patch(Rectangle((j+.08,i-.44),.72,.88,color=co,lw=0))
tipy={t.name:i for i,t in enumerate(tips)}
for n,label,co in side:
 lo=min(tipy[k] for k in leaves[id(n)]);hi=max(tipy[k] for k in leaves[id(n)]);mid=(lo+hi)/2
 br.plot([.05,0,0,.05],[lo-.35,lo-.35,hi+.35,hi+.35],c=co,lw=2)
 br.text(.075,mid,label,fontsize=12.4 if '\n' not in label else 11.7,fontstyle='italic' if '\n' not in label else 'normal',weight='bold',va='center',linespacing=1.5)
for a in [ax,lab,tr,br]:a.set_ylim(65.65,-1.7);a.set_yticks([])
ax.set_xlim(0,mx*1.04);ax.spines[['top','right','left']].set_visible(False);ax.set_xlabel('Substitutions per site',fontsize=11)
lab.set_xlim(0,1);lab.axis('off');tr.set_xlim(0,3);tr.axis('off');br.set_xlim(0,1);br.axis('off')
for j,title in enumerate(['Species','3-gene support','Core split']):tr.text(j+.44,1.006,title,transform=tr.get_xaxis_transform(),rotation=90,va='bottom',ha='center',fontsize=7)
fig.text(.035,.98,'Core-genome phylogeny | type-derived species overlay',fontsize=23,weight='bold',va='top');fig.text(.035,.952,'66 genomes | names transferred from matched BenA + CaM + RPB2 entries',fontsize=11.5,va='top')
output(fig,'core66_type_overlay_corrected')
fig=plt.figure(figsize=(15,10));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off');ax.text(.04,.94,'Core-genome overlay | colours and interpretation',fontsize=23,weight='bold')
for i,sp in enumerate(sorted({ann[k]['assigned_species'] for k in names})):
 y=.85-i*.055;ax.add_patch(Rectangle((.04,y-.01),.018,.02,color=palette[sp]));ax.text(.072,y,sp,fontsize=12,va='center')
for i,(s,co) in enumerate([('Three-gene support: SH-aLRT >=80 and UFBoot >=95','#244f67'),('Three-gene support below either threshold','#d39a35'),('Core group monophyletic, or a single tip','#4b777c'),('Transferred group not monophyletic in the core tree','#b76b40')]):
 y=.85-i*.065;ax.add_patch(Rectangle((.39,y-.01),.018,.02,color=co));ax.text(.422,y,s,fontsize=11,va='center')
ax.scatter(.4,.56,marker='*',s=110,c=GOLD,edgecolor='#674c16');ax.text(.422,.56,'Gold star and [T-G]: documented genome type strain',fontsize=11,va='center')
lines=['Leaf names are inferred from the same genome-derived sample in the three-gene tree.',
 'Large right-hand species labels are shown once, only for monophyletic multi-tip groups in the core tree.',
 'The 40-tip mixed region contains 22 terreus and 18 montoensis three-gene assignments.',
 'Its neutral bracket does not identify a species clade. Singletons retain their leaf labels without a clade claim.',
 'No species-name asterisks are used here: weaker three-gene support is encoded in its own column.',
 'Columns, left to right: type-derived species / three-gene clade support / core monophyly.',
 'Current-tree discordance does not establish biological polyphyly; provenance and orthology gates remain open.',
 'Tree topology, branch lengths and all 66 leaf assignments are unchanged from the preceding release.']
for i,s in enumerate(lines):ax.text(.04,.32-i*.032,s,fontsize=11,va='top')
output(fig,'core66_corrected_legend')
import fitz,subprocess
poppler=Path(r'C:\Users\q287js\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
checks=[]
for stem in ['core66_type_overlay_corrected','core66_corrected_legend']:
 pdf=OUT/'figures'/(stem+'.pdf');d=fitz.open(pdf);assert len(d)==1;p=d[0]
 for b in p.get_text('dict')['blocks']:
  for line in b.get('lines',[]):
   for span in line['spans']:
    x0,y0,x1,y1=span['bbox'];assert x0>=0 and y0>=0 and x1<=p.rect.width+.1 and y1<=p.rect.height+.1,(stem,span['text'])
 subprocess.run([str(poppler),'-r','65','-singlefile','-png',str(pdf),str(OUT/'qa'/stem)],check=True,capture_output=True)
 pix=p.get_pixmap(dpi=300,alpha=False);pix.save(str(OUT/'figures'/(stem+'.png')));del pix;d.close();gc.collect()
 checks.append({'stem':stem,'pages':1,'png_dpi':300,'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()})
save(OUT/'qa/checks.json',checks);shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
(OUT/'README.txt').write_text('CORE LABEL CORRECTION | 20 November 2026\n\nThis supersedes the linear core figure from 2026-11-20_v1_type_monophyly.\nThe repeated large species labels were generated for local pure subgroups, not independent type-qualified clades. The corrected figure retains the same leaf assignments, but labels each globally monophyletic multi-tip core group once. The 40-tip terreus/montoensis region is marked Mixed three-gene assignments. Species-name asterisks are removed; support remains in a dedicated column. Gold stars still mark type genomes.\n\nNo inference or biological reclassification was rerun. Previous three-gene and circular figures remain available in the v1 folder. The mapping is still conditional on the open provenance, locus orthology and core-tree validation gates.\n',encoding='utf-8')
print(json.dumps(record,indent=2))
