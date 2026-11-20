"""Preserve original taxon labels and explicitly restore all type-anchor information."""
from pathlib import Path
import sys,json,hashlib,collections,gc,shutil,os,re
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
limits=constrain_process()
from terreus_tree_helpers import parse_newick
OLD=ROOT/'stakeholder_release/2026-11-20_v1_type_monophyly'
OUT=ROOT/'stakeholder_release/2026-11-20_v3_original_taxon_labels'
for d in ['figures','data','qa','scripts']:(OUT/d).mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(OLD/'qa/mpl')
def load(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
treepath=OLD/'data/core66_display_rooted.treefile';r=parse_newick(treepath.read_text())
audit=load(OLD/'data/assignment_audit.json');ann={x['tip']:x for x in audit}
original_ann=load(ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene/data/new_tree_annotations.json')
source_metadata_path=ROOT/'review_reports/2026-11-10/ncbi_core_tree_metadata/core_tree_metadata_complete.json'
source_metadata={x['summary']['sample_id']:x['summary'] for x in load(source_metadata_path)}
for a in audit:
 if not a['in_core66']:continue
 a['previous_plot_taxon']=a['reported_species']
 organism=source_metadata[a['sample']].get('ncbi_organism','')
 a['original_source_organism']=organism
 if organism:
  m=re.match(r'^Aspergillus ([a-z][a-z-]+)(?:\s|$)',organism)
  assert m,(a['sample'],organism)
  a['reported_species']='A. '+m.group(1)
  a['original_taxon_basis']='Public organism field in archived NCBI source metadata; species binomial only, strain suffix retained separately.'
 else:a['original_taxon_basis']='Original project taxon label; no public organism record established for this study sample.'
 a['reported_name_changed']=a['assigned_species']!=a['reported_species']
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
record['correction']='Each tip shows inferred species (original taxon) | original strain ID. All three actual genome types retain gold stars and their full documented culture identifiers. All nine taxon type anchors are explicitly listed in the shared legend; references without core genomes are not inserted into the core topology.'
record['label_format']='Inferred species (original reported taxon) | original strain ID | TYPE culture, where applicable'
record['tip_labels_unchanged']=False
record['species_assignments_unchanged']=True
record['original_taxon_source']={'path':str(source_metadata_path),'sha256':hashlib.sha256(source_metadata_path.read_bytes()).hexdigest()}
record['earlier_plot_taxa_corrected_to_source']=[{'sample':a['sample'],'earlier_plot_taxon':a['previous_plot_taxon'],'original_source_taxon':a['reported_species']} for a in audit if a['in_core66'] and a['reported_species']!=a['previous_plot_taxon']]
record['original_taxon_counts']=dict(collections.Counter(a['reported_species'] for a in audit if a['in_core66']))
def tip_label(k):
 a=ann[k];o=original_ann[k]
 s=a['assigned_species']+' ('+a['reported_species']+') | '+o['strain']
 if a['is_type']:s+=' | TYPE: '+o['genome_type_culture']
 return s
for a in audit:
 if a['in_core66']:a['display_label']=tip_label(a['tip'])
save(OUT/'data/correction_audit.json',record);save(OUT/'data/core_tip_assignments.json',[a for a in audit if a['in_core66']])
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
INK='#203849';GRAY='#70818f';GOLD='#d9a72f';SUP={'supported':'#244f67','topology_only':'#d39a35'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','text.color':INK})
def output(fig,stem):
 fig.savefig(OUT/'figures'/(stem+'.pdf'),facecolor='white');fig.savefig(OUT/'figures'/(stem+'.svg'),facecolor='white');plt.close(fig);gc.collect()
fig=plt.figure(figsize=(24,16));ax=fig.add_axes([.03,.055,.265,.88]);lab=fig.add_axes([.31,.055,.395,.88],sharey=ax);tr=fig.add_axes([.72,.055,.045,.88],sharey=ax);br=fig.add_axes([.80,.055,.18,.88],sharey=ax);mx=max(xs.values())
def branches(n):
 if n.children:
  ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],c=GRAY,lw=.65)
  for c in n.children:
   ss={ann[k]['assigned_species'] for k in leaves[id(c)]};co=palette[next(iter(ss))] if len(ss)==1 else GRAY
   ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,c=co,lw=.85);branches(c)
 else:ax.plot([xs[id(n)],mx*1.025],[ys[id(n)]]*2,c='#dce1e5',lw=.4,ls=':')
branches(r)
for i,t in enumerate(tips):
 a=ann[t.name];label=tip_label(t.name)
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
fig.text(.03,.98,'Core-genome phylogeny | inferred species and original taxon',fontsize=23,weight='bold',va='top');fig.text(.03,.952,'66 genomes | inferred species (original taxon) | strain identifier | documented types explicitly marked',fontsize=11.5,va='top')
output(fig,'core66_inferred_and_original_linear')
# Same names and type identifiers in the circular view; topology shown as a cladogram.
import numpy as np,math
depth={id(r):0}
def depths(n):
 for child in n.children:depth[id(child)]=depth[id(n)]+1;depths(child)
depths(r);md=max(depth.values());step=2*math.pi/66
fig,ax=plt.subplots(figsize=(27,27),subplot_kw={'projection':'polar'})
def circle(n):
 radius=depth[id(n)]/md;angle=ys[id(n)]*step
 if n.children:
  aa=[ys[id(c)]*step for c in n.children];ax.plot(np.linspace(min(aa),max(aa),50),[radius]*50,c=GRAY,lw=.65)
  for child in n.children:ax.plot([ys[id(child)]*step]*2,[radius,depth[id(child)]/md],c=GRAY,lw=.65);circle(child)
 else:
  a=ann[n.name];ax.plot([angle]*2,[radius,1.015],c='#cbd4dc',lw=.4,ls=':')
  cols=[palette[a['assigned_species']],SUP[a['confidence']],'#b76b40' if a['core_monophyly'] is False else '#4b777c']
  for j,co in enumerate(cols):ax.bar(angle,.031,width=step*.94,bottom=1.035+j*.045,color=co,lw=0)
  if a['is_type']:ax.scatter(angle,1.19,marker='*',s=65,c=GOLD,edgecolor='#674c16',lw=.3)
  deg=math.degrees(angle);left=90<deg<270
  ax.text(angle,1.225,tip_label(n.name),rotation=deg+180 if left else deg,rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=9,weight='bold' if a['is_type'] else 'normal')
circle(r);ax.set_ylim(0,2.45);ax.axis('off');fig.subplots_adjust(left=.015,right=.985,bottom=.02,top=.93)
fig.text(.03,.98,'Core-genome phylogeny | inferred species and original taxon',fontsize=25,weight='bold',va='top');fig.text(.03,.958,'66 genomes | circular cladogram | parentheses retain original taxon labels',fontsize=12,va='top')
output(fig,'core66_inferred_and_original_circular')
fig=plt.figure(figsize=(18,13));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off');ax.text(.04,.95,'Species, original labels and type anchors',fontsize=23,weight='bold')
ax.text(.04,.905,'Inferred taxon and its three-gene type anchor',fontsize=12,weight='bold')
for i,sp in enumerate(sorted({ann[k]['assigned_species'] for k in names})):
 y=.86-i*.052;ax.add_patch(Rectangle((.04,y-.009),.014,.018,color=palette[sp]));ax.text(.065,y,sp,fontsize=11.5,va='center')
 g=next(g for g in groups if g['species']==sp)
 reference_anchor=next(k for k in g['anchors'] if not k.startswith('LOCAL_'))
 culture=original_ann[reference_anchor]['strain']
 actual=next((original_ann[k].get('genome_type_culture') for k in g['anchors'] if k in names and original_ann[k].get('genome_type_culture')),None)
 ax.text(.25,y,culture+(' [also a core type genome]' if actual else ' [three-gene reference only]'),fontsize=10.5,va='center')
for i,(s,co) in enumerate([('Three-gene support: SH-aLRT >=80 and UFBoot >=95','#244f67'),('Three-gene support below either threshold','#d39a35'),('Core group monophyletic, or a single tip','#4b777c'),('Transferred group not monophyletic in the core tree','#b76b40')]):
 y=.86-i*.06;ax.add_patch(Rectangle((.63,y-.009),.014,.018,color=co));ax.text(.655,y,s,fontsize=9.5,va='center')
ax.scatter(.64,.59,marker='*',s=110,c=GOLD,edgecolor='#674c16');ax.text(.66,.59,'Gold star + TYPE culture: actual core type genome',fontsize=10,va='center')
ax.text(.63,.535,'Example label',fontsize=12,weight='bold')
ax.text(.63,.505,'A. montoensis (A. terreus) | ATCC20541',fontsize=10.5)
ax.text(.63,.475,'Inferred species (original taxon) | original strain ID',fontsize=10)
lines=['Names outside parentheses are inferred from the same genome-derived sample in the three-gene tree.',
 'Parentheses use the public source organism taxon; study isolates without a public record retain their original project taxon.',
 'Actual core type genomes: SRR10507650 = CBS 116.37; SRR35107300 = IBT 26384 = CBS 124230;',
 'SRR8399223 = CBS 503.65. Their leaves show gold stars, bold original strain IDs and full type culture identifiers.',
 'Other listed types anchor assignments in the three-gene tree only; no artificial type tips have been added to the core tree.',
 'Large right-hand species labels are shown once, only for monophyletic multi-tip groups in the core tree.',
 'The 40-tip mixed region contains 22 terreus and 18 montoensis three-gene assignments.',
 'Its neutral bracket does not identify a species clade. Singletons retain their leaf labels without a clade claim.',
 'No species-name asterisks are used here: weaker three-gene support is encoded in its own column.',
 'Linear columns and circular rings (inner to outer): type-derived species / three-gene support / core monophyly.',
 'Current-tree discordance does not establish biological polyphyly; provenance and orthology gates remain open.',
 'All 66 assignments and the core topology are unchanged. Linear branches retain lengths; the circular view is a cladogram.']
for i,s in enumerate(lines):ax.text(.04,.355-i*.025,s,fontsize=10.5,va='top')
output(fig,'core66_original_labels_and_types_legend')
import fitz,subprocess
poppler=Path(r'C:\Users\q287js\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
checks=[]
for stem in ['core66_inferred_and_original_linear','core66_inferred_and_original_circular','core66_original_labels_and_types_legend']:
 pdf=OUT/'figures'/(stem+'.pdf');d=fitz.open(pdf);assert len(d)==1;p=d[0]
 for b in p.get_text('dict')['blocks']:
  for line in b.get('lines',[]):
   for span in line['spans']:
    x0,y0,x1,y1=span['bbox'];assert x0>=0 and y0>=0 and x1<=p.rect.width+.1 and y1<=p.rect.height+.1,(stem,span['text'])
 subprocess.run([str(poppler),'-r','65','-singlefile','-png',str(pdf),str(OUT/'qa'/stem)],check=True,capture_output=True)
 pix=p.get_pixmap(dpi=300,alpha=False);pix.save(str(OUT/'figures'/(stem+'.png')));del pix;d.close();gc.collect()
 checks.append({'stem':stem,'pages':1,'png_dpi':300,'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()})
save(OUT/'qa/checks.json',checks);shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
(OUT/'README.txt').write_text('ORIGINAL TAXON LABELS AND TYPE ANCHORS | 20 November 2026\n\nCore leaf format: inferred species (original reported taxon) | original strain ID. Original taxon is taken from the public source organism field where available; study isolates without a public organism record retain their original project taxon. Earlier plot names are preserved separately in the audit. Five formerly displayed pseudoterreus labels (AT2022, IFO6365, IFO6365_DRR204676, TN-484, TN-484_DRR191518) have terreus in the original public source record; their parentheses now correctly show terreus.\n\nThe three documented type-genome leaves show TYPE and their full culture identifiers, in bold with gold stars. The single shared legend lists all nine type anchors, distinguishes three-gene-only type references, and explains every colour.\n\nIncludes linear and circular single-page vector PDFs, SVGs and 300-dpi PNGs. No artificial type-reference tip has been inserted into the core topology. All 66 inferred species assignments and the underlying topology are unchanged. The three-gene tree files from v1 are preserved.\n',encoding='utf-8')
print(json.dumps(record,indent=2))
