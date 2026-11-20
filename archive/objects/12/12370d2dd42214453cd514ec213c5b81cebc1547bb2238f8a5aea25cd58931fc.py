from pathlib import Path
import sys, os, json, math, shutil, hashlib, textwrap, gc
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
os.environ['MPLCONFIGDIR']=str(OUT/'qa'/'matplotlib_cache')
sys.path.insert(0,str(ROOT/'comparative_visualization/2026-11-13/vendor'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pymupdf
from terreus_tree_helpers import parse_newick,Node
P=json.loads((OUT/'data/package.json').read_text());F=OUT/'figures';Q=OUT/'qa';OLD=ROOT/'comparative_visualization/2026-11-11';NOV13=ROOT/'comparative_visualization/2026-11-13'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':'#213c50','axes.labelcolor':'#213c50','svg.fonttype':'none','pdf.fonttype':42})
source_hashes={}
def source(p):
 p=Path(p);source_hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return p
def save(fig,stem,pdf=None):
 if pdf:pdf.savefig(fig,bbox_inches='tight',pad_inches=.3)
 else:fig.savefig(F/(stem+'.pdf'),bbox_inches='tight',pad_inches=.3)
 for ext in ['svg','png']:fig.savefig(F/(stem+'.'+ext),dpi=100,bbox_inches='tight',pad_inches=.3)
 plt.close(fig);gc.collect()
def page(title,subtitle=''):
 fig=plt.figure(figsize=(16.54,11.69));fig.text(.04,.95,title,fontsize=22,weight='bold');fig.text(.04,.91,subtitle,fontsize=11,color='#526778');return fig
def swatch(fig,x,y,color,label,size=11):
 fig.add_artist(Rectangle((x,y-.006),.012,.014,transform=fig.transFigure,facecolor=color,edgecolor='#88929a',lw=.4));fig.text(x+.019,y,label,va='center',fontsize=size)
palette=P['palette']

# Preserve the latest type-marked topology/layout exactly, and attach complete colour keys.
marker_source=source(NOV13/'marker_provenance/combined_type_marked/kombinierter_baum_typstaemme.pdf')
with PdfPages(Q/'marker_keys.pdf') as pdf:
 fig=page('Three-gene tree | colour and symbol key','Applies to both views: circular rings from inside out; vertical tracks from left to right.')
 lines=[
 '1  Source: teal = genome-extracted [G]; orange = taxonomic reference [R] or historical 2011 reference [R11].',
 '2  Loci: number of available markers in the BenA + CaM + RPB2 panel; missing sequence remains missing.',
 '3  Species: reviewed taxon; grey = unresolved, lavender-grey = source-name conflict. Labels retain reported names.',
 '4  Reference type clade: tree-specific, containing an ex-type and SH-aLRT >=80 plus UFBoot >=95.',
 '5  Genome type: documented ex-type genome and supported compatible in-tree group; blank = no such supported clade.',
 'Gold star + [TYP] = reference ex-type entry; gold star + [TYP-G] = documented ex-type genome.',
 'Black internal-node dot = both support thresholds passed. * after a taxon = unresolved reported name.',
 '283 tips = 76 genome-derived + 207 reference entries. 31 reference type entries + 3 genome type anchors are not 34 unique cultures.'
 ]
 for i,t in enumerate(lines):fig.text(.04,.865-i*.029,t,fontsize=11)
 for x,col,label in [(.04,'#087f8c','Genome source'),(.25,'#bd642d','Reference source'),(.47,'#d7794c','1 locus'),(.62,'#e6ad51','2 loci'),(.77,'#408f87','3 loci')]:swatch(fig,x,.605,col,label)
 fig.text(.04,.565,'Species colours are shared across the species and type-clade tracks',fontsize=14,weight='bold')
 for i,(tax,col) in enumerate(palette.items()):swatch(fig,.04+(i//18)*.49,.531-(i%18)*.024,col,tax)
 fig.text(.04,.06,'Grey mixed-type clades contain more than one type species; they are not species assignments. White does not mean absence.',fontsize=10)
 fig.text(.04,.033,'Circular view is a cladogram. The vertical view preserves uncompressed branch lengths. No unverified sensu-stricto or subspecies labels added.',fontsize=10)
 save(fig,'01_three_gene_colour_key',pdf)
 fig=page('Three-gene tree | reference type-clade codes','Codes are specific to this tree. Nested or mixed-type clades are retained without forcing a species call.')
 fig.text(.04,.86,'Code',weight='bold');fig.text(.11,.86,'Type species within clade',weight='bold');fig.text(.86,.86,'SH-aLRT / UFBoot',weight='bold')
 y=.825
 for c in P['type_clades']:
  tax=c['type_species'];col=palette.get(tax,palette['Mixed type species']);swatch(fig,.04,y,col,c['clade_id'],10)
  wrapped=textwrap.wrap(tax,112);fig.text(.11,y,'\n'.join(wrapped),fontsize=10,va='center');fig.text(.88,y,f"{c['SH_aLRT']:g} / {c['UFBoot']:g}",fontsize=10,va='center')
  y-=.024+(.011 if len(wrapped)>1 else 0)
 fig.text(.04,.10,'Documented genome types: A. floccosus CBS 116.37 (SRR10507650); A. aureoterreus CBS 503.65 (SRR8399223);',fontsize=10)
 fig.text(.04,.075,'A. hortae IBT 26384 = CBS 124230 (SRR35107300). Six conflicting BioSample type claims remain withheld.',fontsize=10)
 fig.text(.04,.045,'The full per-tip annotations and ex-type entry list are included in data/. Source provenance and type status are separate variables.',fontsize=10)
 save(fig,'01_three_gene_type_clade_key',pdf)
doc=pymupdf.open(marker_source);keys=pymupdf.open(Q/'marker_keys.pdf');doc.insert_pdf(keys);doc.set_metadata({'title':'Aspergillus terreus group | BenA + CaM + RPB2','subject':'Latest combined tree with type strains, sequence provenance and full colour-coded keys'});doc.save(F/'01_three_gene_tree.pdf');doc.close();keys.close()
for suffix,new in [('kreisfoermig','circular'),('linear','vertical')]:
 for ext in ['svg','png']:shutil.copy2(source(marker_source.parent/f'kombinierter_baum_typstaemme_{suffix}.{ext}'),F/f'01_three_gene_{new}.{ext}')
shutil.copy2(source(marker_source.parent/'markierte_typstaemme.tsv'),OUT/'data/type_strain_entries.tsv')

# Exact 66-tip core topology, displayed with the same outgroup split as the previous release.
root=parse_newick(source(OLD/'sources/input_trees/core_genome_original.treefile').read_text())
for t in root.leaves():t.name='LOCAL_'+t.name.removesuffix('_transcripts').removesuffix('.transcripts')
keep={t.name for t in root.leaves()}
def splits(n):
 out=set()
 def walk(x):
  a={x.name} if not x.children else set().union(*(walk(c) for c in x.children));b=keep-a
  if len(a)>1 and len(b)>1:out.add(min((tuple(sorted(a)),tuple(sorted(b))),key=lambda z:(len(z),z)))
  return a
 walk(n);return out
before=splits(root)
def reroot(root,name):
 graph={};nodes={}
 def walk(n):
  nodes[id(n)]=n;graph.setdefault(id(n),[])
  for c in n.children:graph[id(n)].append((id(c),c.length));graph.setdefault(id(c),[]).append((id(n),c.length));walk(c)
 walk(root);tip=next(id(n) for n in root.leaves() if n.name==name);nei,dist=graph[tip][0]
 def branch(i,parent,length):
  kids=[branch(j,i,d) for j,d in graph[i] if j!=parent];n=Node(nodes[i].name if not kids else '',length,kids)
  if len(kids)==1:kids[0].length+=length;return kids[0]
  return n
 return Node('',0,[Node(name,dist/2),branch(nei,tip,dist/2)])
root=reroot(root,'LOCAL_SRR10507650');assert splits(root)==before
ann={a['tip']:a for a in P['core_annotations']};metrics={'LOCAL_'+m['sample']:m for m in P['genomes']};assert keep==set(ann)
cohortcols={'Environmental study':'#087f8c','Project R22-R24':'#bc7fa2','Other genome':'#bdc7d0'}
mgcols={'MG1':'#3b69b1','MG2':'#bd7b30','MG3':'#7853a5','unknown':'#eeeeee'}
lovcols={'nine_gene_block':'#347f73','partial_or_divergent_block':'#d5a64b','no_convincing_colocalized_block':'#b4bac1'}
norm=LogNorm(1e3,1e7);cmap=plt.get_cmap('Blues')
def tracks(k):
 a=ann[k];m=metrics[k]
 return [cohortcols[m['project_group']],mgcols[m['major_genotype']],palette[a['species_ring']],palette[a['type_clade_category']],palette.get(a.get('genome_type_species',''),'#f3f4f6'),lovcols[m['block_status']],cmap(norm(m['N50_bp']))]
def label(k):return ann[k]['display_label'].replace('[TG]','[TYP-G]')
def positions(clado):
 tips=root.leaves();y={id(t):i for i,t in enumerate(tips)};x={id(root):0}
 def walk(n):
  for c in n.children:x[id(c)]=x[id(n)]+(1 if clado else c.length);walk(c)
  if n.children:y[id(n)]=sum(y[id(c)] for c in n.children)/len(n.children)
 walk(root);return x,y,tips
with PdfPages(F/'02_core_genome_tree.pdf') as pdf:
 fig,ax=plt.subplots(figsize=(26,25),subplot_kw={'projection':'polar'});x,y,tips=positions(True);mx=max(x.values());step=2*math.pi/len(tips)
 def draw(n):
  angle=y[id(n)]*step;r=x[id(n)]/mx
  if n.children:
   aa=[y[id(c)]*step for c in n.children];ax.plot(np.linspace(min(aa),max(aa),100),[r]*100,color='#778592',lw=.8)
   for c in n.children:ax.plot([y[id(c)]*step]*2,[r,x[id(c)]/mx],color='#778592',lw=.8);draw(c)
  else:
   k=n.name;m=metrics[k];ax.plot([angle]*2,[r,1.025],color='#cbd3da',ls=':',lw=.6)
   ax.scatter(angle,1.025,s=20,color=cohortcols[m['project_group']],edgecolor='#526778',linewidth=.4,zorder=4)
   for j,col in enumerate(tracks(k)):ax.bar(angle,.026,width=step*.94,bottom=1.055+j*.036,color=col,lw=0)
   if ann[k].get('genome_type_species'):ax.scatter(angle,1.324,marker='*',s=85,color='#e8af27',edgecolor='#533f12',linewidth=.5,zorder=6)
   deg=math.degrees(angle);left=90<deg<270;ax.text(angle,1.35,label(k),rotation=deg+(180 if left else 0),rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=11,color='#087f8c' if m['project_group']=='Environmental study' else '#6e3f64' if m['project_group']=='Project R22-R24' else '#34495b',weight='bold' if m['project_group']!='Other genome' else 'normal')
 draw(root);ax.set_ylim(0,1.99);ax.axis('off');fig.subplots_adjust(left=.01,right=.99,bottom=.075,top=.91)
 fig.text(.035,.972,'Aspergillus terreus group | core-genome context',fontsize=27,weight='bold')
 fig.text(.035,.945,'66 genomes | 4,661 single-copy loci | original GTR+F+R10 topology | EXPLORATORY: no bootstrap or SH-aLRT',fontsize=14)
 fig.text(.035,.925,'Klebsiella-style circular layout, study markers and metadata rings. Cladogram: radial distances do not represent branch lengths.',fontsize=12)
 fig.text(.035,.045,'Rings, inside out: 1 cohort | 2 major genotype | 3 reviewed species | 4 projected reference type clade | 5 genome type anchor | 6 lovastatin homology | 7 N50',fontsize=12)
 fig.text(.035,.025,'Gold star = documented ex-type genome. * = unresolved reported species. Projection from the marker tree is not support for a core-genome clade. Full colour key follows.',fontsize=11)
 save(fig,'02_core_genome_circular',pdf)
 x,y,tips=positions(False);mx=max(x.values());fig=plt.figure(figsize=(25,20));ax=fig.add_axes([.04,.09,.38,.80]);lab=fig.add_axes([.435,.09,.355,.80],sharey=ax);tr=fig.add_axes([.80,.09,.18,.80],sharey=ax)
 def drawlin(n):
  if n.children:
   ax.plot([x[id(n)]]*2,[min(y[id(c)] for c in n.children),max(y[id(c)] for c in n.children)],c='#778592',lw=.8)
   for c in n.children:ax.plot([x[id(n)],x[id(c)]],[y[id(c)]]*2,c='#778592',lw=.8);drawlin(c)
  else:ax.plot([x[id(n)],mx*1.01],[y[id(n)]]*2,c='#cbd3da',ls=':',lw=.6)
 drawlin(root)
 for i,t in enumerate(tips):
  m=metrics[t.name];lab.text(0,i,label(t.name),fontsize=11,va='center',color='#087f8c' if m['project_group']=='Environmental study' else '#6e3f64' if m['project_group']=='Project R22-R24' else '#34495b',weight='bold' if m['project_group']!='Other genome' else 'normal')
  if ann[t.name].get('genome_type_species'):lab.scatter(-.035,i,marker='*',s=70,color='#e8af27',edgecolor='#533f12',linewidth=.5,clip_on=False)
  for j,col in enumerate(tracks(t.name)):tr.add_patch(Rectangle((j-.4,i-.44),.8,.88,facecolor=col,lw=0))
 for a in [ax,lab,tr]:a.set_ylim(len(tips)-.3,-1.5);a.set_yticks([])
 ax.set_xlim(0,mx*1.02);ax.spines[['top','left','right']].set_visible(False);ax.set_xlabel('Substitutions per site | uncompressed')
 lab.set_xlim(0,1);lab.axis('off');tr.set_xlim(-.6,6.6);tr.set_xticks(range(7),['1 Cohort','2 MG','3 Species','4 Ref. clade','5 Type anchor','6 Lov block','7 N50'],rotation=40,ha='left');tr.xaxis.tick_top();tr.spines[:].set_visible(False)
 fig.text(.04,.965,'Core-genome tree | uncompressed branch lengths',fontsize=25,weight='bold')
 fig.text(.04,.94,'Original 66-tip topology; display rooting on SRR10507650 preserves all unrooted splits. Rooting remains a scientific gate.',fontsize=12)
 fig.text(.04,.92,'EXPLORATORY: no branch-support estimates. The same seven metadata tracks are used in the circular view.',fontsize=12,color='#a36c23')
 fig.text(.04,.045,'* Reported species unresolved by the three markers. Gold star / [TYP-G]: documented ex-type genome. Colour key on next page.',fontsize=11)
 save(fig,'02_core_genome_vertical',pdf)
 fig=page('Core-genome tree | metadata key','Same fixed species palette as the marker tree. No supported core type-clade assignments are inferred from this topology.')
 fig.text(.04,.85,'1  Cohort',fontsize=14,weight='bold')
 for i,(label_,col) in enumerate(cohortcols.items()):swatch(fig,.04,.81-i*.035,col,label_+' (metadata-based)')
 fig.text(.04,.665,'2  Major genotype',fontsize=14,weight='bold')
 for i,(label_,col) in enumerate(mgcols.items()):swatch(fig,.04,.625-i*.035,col,label_.replace('unknown','Unknown / not recorded'))
 fig.text(.04,.44,'6  Lovastatin homology state',fontsize=14,weight='bold')
 for i,(key,txt) in enumerate([('nine_gene_block','Nine-gene co-localised block'),('partial_or_divergent_block','Partial / divergent block'),('no_convincing_colocalized_block','No convincing co-localised block')]):swatch(fig,.04,.40-i*.035,lovcols[key],txt)
 fig.text(.04,.27,'These states do not demonstrate production or true absence.',fontsize=10)
 fig.text(.04,.225,'7  Assembly N50 (bp; log scale)',fontsize=14,weight='bold')
 cbax=fig.add_axes([.04,.16,.34,.025]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbax,orientation='horizontal',ticks=[1e3,1e4,1e5,1e6,1e7])
 fig.text(.5,.85,'3–5  Reviewed species, projected clade, type anchor',fontsize=14,weight='bold')
 cats=sorted(set(a['species_ring'] for a in ann.values())|set(a['type_clade_category'] for a in ann.values())|set(a.get('genome_type_species','') for a in ann.values())- {''})
 for i,tax in enumerate(cats):swatch(fig,.5,.805-i*.032,palette.get(tax,'#f3f4f6'),tax)
 y=.805-len(cats)*.032-.025;swatch(fig,.5,y,'#f3f4f6','No documented genome-type anchor at this tip')
 fig.text(.5,y-.09,'Track 4 is projected from the combined reference tree.\nTrack 5 marks the three documented ex-type genomes only.\nNeither track establishes supported clades in the core tree.',fontsize=11,linespacing=1.6)
 fig.text(.04,.065,'Project metadata: 12 environmental isolates plus R22/R23/R24. Other genome does not imply complete provenance; missing MG remains unknown.',fontsize=10)
 fig.text(.04,.035,'A high N50 is not completeness. Raw-read QC, harmonised BUSCO, rooting, support and annotation/model sensitivity remain separate gates.',fontsize=10)
 save(fig,'02_core_genome_colour_key',pdf)

# Preserve exact latest architecture and genome-wide composition figures.
for old,new in [('lovastatin_plain_linear_architecture','03_lovastatin_linear')]:
 for ext in ['pdf','svg','png']:shutil.copy2(source(NOV13/'plain_views'/f'{old}.{ext}'),F/f'{new}.{ext}')
shutil.copy2(source(NOV13/'plain_views/region_features.tsv'),OUT/'data/lovastatin_region_features.tsv')

# Reuse the validated Circos layout with the tree taxon palette and the shared N50 scale.
# Extract only the two plotting functions, avoiding the original module's output/cache setup.
import ast
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from pycirclize import Circos
layout=source(ROOT/'plot_plain_terreus_views.py')
nodes=[n for n in ast.parse(layout.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in ['bar','genome_view']]
order=json.loads(source(NOV13/'data/figure_order.json').read_text());ga=json.loads(source(NOV13/'data/sample_annotations.json').read_text());gm=json.loads(source(NOV13/'data/genome_metrics.json').read_text())
by_sample={s:m for m in gm.values() for s in m['samples']};taxa=sorted({ga['LOCAL_'+s]['taxon'] for s in order})
def genome_label(s):
 a=ga['LOCAL_'+s];flag='*' if a['species_ring'] in ('Unresolved','Name conflict') else '';return f'{order.index(s)+1:02d}  {a["taxon"]}{flag} | {s}'
env={'plt':plt,'math':math,'Circos':Circos,'Normalize':Normalize,'LogNorm':LogNorm,'ScalarMappable':ScalarMappable,'Patch':Patch,'order':order,'ann':ga,'by_sample':by_sample,'taxa':taxa,'taxcolors':{t:palette.get(t,'#8d969f') for t in taxa},'label':genome_label,'save':lambda fig,stem:save(fig,'04_whole_genome_circos')}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(layout),'exec'),env);env['genome_view']()
for fn in ['genome_metrics.json','figure_order.json','sample_annotations.json']:shutil.copy2(source(NOV13/'data'/fn),OUT/'data'/fn)

# Re-present the existing, matched-cohort statistical results without rerunning inference.
fig,(ax,bx)=plt.subplots(1,2,figsize=(17,8),gridspec_kw={'width_ratios':[1.12,1]});fig.subplots_adjust(left=.17,right=.96,bottom=.19,top=.77,wspace=.30)
boot={a['marker_choice']:a for a in P['bootstrap']};cols={r['marker_choice']:plt.get_cmap('tab10')(i) for i,r in enumerate(P['ranking'])};cols['Core genome']='#253e59'
for i,r in enumerate(P['ranking']):
 b=boot[r['marker_choice']];col=cols[r['marker_choice']];ax.plot([b['percentile_2_5'],b['percentile_97_5']],[i,i],color=col,lw=3,alpha=.65);ax.scatter(r['normalized_RF_to_core'],i,s=70,color=col,zorder=3);ax.text(.898,i,f"{r['normalized_RF_to_core']:.3f}",fontsize=10,va='center')
ax.set_yticks(range(7),[r['marker_choice'].replace('_',' + ') for r in P['ranking']]);ax.invert_yaxis();ax.set_xlim(.57,.94);ax.set_xlabel('Normalised Robinson–Foulds distance to core tree');ax.set_title('A  Marker-tree proximity',loc='left',weight='bold');ax.grid(axis='x',alpha=.16);ax.spines[['top','right']].set_visible(False)
offsets={'Core genome':(9,-16),'BenA':(-8,12),'CaM':(8,-16),'RPB2':(-50,10),'BenA_CaM':(-45,-18),'BenA_RPB2':(-48,14),'CaM_RPB2':(8,7),'BenA_CaM_RPB2':(-60,-22)}
nm=[a for a in P['ordination'] if a['method']=='NMDS']
for r in nm:
 bx.scatter(r['axis1'],r['axis2'],marker='*' if r['tree']=='Core genome' else 'o',s=170 if r['tree']=='Core genome' else 70,color=cols[r['tree']]);bx.annotate(r['tree'].replace('_',' + '),(r['axis1'],r['axis2']),xytext=offsets[r['tree']],textcoords='offset points',fontsize=10)
bx.set(xlim=(-1.05,1.24),ylim=(-.62,.80),xlabel='NMDS axis 1',ylabel='NMDS axis 2');bx.set_aspect('equal');bx.set_title(f"B  Eight-tree NMDS | stress-1 = {nm[0]['stress1']:.4f}",loc='left',weight='bold');bx.spines[['top','right']].set_visible(False);bx.axhline(0,color='#e4e8eb',lw=.7);bx.axvline(0,color='#e4e8eb',lw=.7)
fig.text(.035,.945,'Which marker combination most resembles the core tree?',fontsize=22,weight='bold');fig.text(.035,.90,'All comparisons use the same 66 tips. BenA + CaM + RPB2 is closest overall (nRF 0.635), with substantial remaining discordance.',fontsize=12)
fig.text(.035,.105,'Dots: ML-tree distance. Bars: 2.5th–97.5th percentiles of 1,000 marker bootstrap trees per choice versus the fixed core tree.',fontsize=11)
fig.text(.035,.065,'The core tree has no bootstrap / SH-aLRT analysis. Similarity is not validated taxonomic accuracy; overlapping intervals are not pairwise significance tests.',fontsize=10)
fig.text(.035,.03,'NMDS summarises pairwise unweighted RF distances among eight trees. The 2D view is approximate; the complete distance matrix is supplied.',fontsize=10)
save(fig,'05_marker_comparison')
shutil.copy2(source(OLD/'core66_combinations/normalized_RF_matrix.tsv'),OUT/'data/normalized_RF_matrix.tsv')

with PdfPages(Q/'cover.pdf') as pdf:
 fig=page('Aspergillus terreus group | stakeholder figure set','Release 2026-11-17_v1 | latest reviewed analyses and figure exports | no sequence inference rerun')
 fig.text(.04,.84,'A complete visual package with explicit scientific gates',fontsize=17,weight='bold')
 fig.text(.04,.785,'66 core-tree genomes • 77 assembly labels / 72 canonical sequences • 283 reference-inclusive three-gene tips',fontsize=12)
 y=.71
 for ident,title,fn,scope,purpose,caveat in P['figures']:
  fig.text(.04,y,ident+'  '+title,fontsize=14,weight='bold');fig.text(.08,y-.027,scope,fontsize=11);fig.text(.08,y-.054,caveat,fontsize=10,color='#756046');y-=.105
 fig.text(.04,.135,'Read with analysis_gates.xlsx: computation completed and scientific gate passed are deliberately separate statuses.',fontsize=11,weight='bold')
 fig.text(.04,.095,'Held back: legacy SNP, secretome and functional clustering figures. Ancestral-state / evo-devo analyses are proposed, not run.',fontsize=11)
 fig.text(.04,.055,'Large-format vector pages preserve all leaf and gene labels. Zoom or print at native size; each main figure is also supplied separately.',fontsize=10)
 pdf.savefig(fig);plt.close(fig)
book=pymupdf.open();book.insert_pdf(pymupdf.open(Q/'cover.pdf'));toc=[[1,'Release overview',1]]
for _,title,fn,*_ in P['figures']:
 toc.append([1,title,len(book)+1]);d=pymupdf.open(F/fn);book.insert_pdf(d);d.close()
book.set_toc(toc);book.set_metadata({'title':'Aspergillus terreus group | stakeholder figure set 2026-11-17_v1'});book.save(OUT/'stakeholder_figures.pdf');book.close()

original=pymupdf.open(marker_source);exported=pymupdf.open(F/'01_three_gene_tree.pdf')
preserved=all(original[i].get_text()==exported[i].get_text() and original[i].get_pixmap(matrix=pymupdf.Matrix(.15,.15)).samples==exported[i].get_pixmap(matrix=pymupdf.Matrix(.15,.15)).samples for i in range(2))
original.close();exported.close()
checks={'core_tip_count':len(keep)==66,'core_splits_preserved':splits(root)==before,'core_metadata_complete':keep==set(ann),'core_metrics_complete':keep<=set(metrics),'all_77_labels_in_metrics':len(P['genomes'])==77,'unique_sequences_72':len({m['sequence_sha256'] for m in P['genomes']})==72,'latest_combined_first_two_pages_preserved':preserved}
for fn in [x[2] for x in P['figures']]:
 doc=pymupdf.open(F/fn)
 for i,p in enumerate(doc):
  scale=min(1.3,1500/p.rect.width,1800/p.rect.height);p.get_pixmap(matrix=pymupdf.Matrix(scale,scale)).save(str(Q/f'{Path(fn).stem}_page{i+1}.png'))
  if fn=='01_three_gene_tree.pdf' and i==1:
   # Inspect exact vertical track region and a representative label section at native resolution.
   p.get_pixmap(matrix=pymupdf.Matrix(1,1),clip=pymupdf.Rect(p.rect.width*.76,p.rect.height*.03,p.rect.width*.99,p.rect.height*.22)).save(str(Q/'marker_vertical_tracks_detail.png'))
 doc.close()
checks['all_source_files_unchanged']=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in source_hashes.items())
assert all(checks.values()),checks
(OUT/'data/figure_source_hashes.json').write_text(json.dumps(source_hashes,indent=2))
(OUT/'qa/figure_validation.json').write_text(json.dumps({'checks':checks,'resource_limits':LIMITS},indent=2))
print('FIGURES COMPLETE',json.dumps(checks),flush=True)
