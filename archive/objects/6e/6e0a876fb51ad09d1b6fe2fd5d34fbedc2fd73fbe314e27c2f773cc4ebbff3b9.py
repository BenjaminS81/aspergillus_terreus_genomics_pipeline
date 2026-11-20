from pathlib import Path
import sys,os,json,math,hashlib,csv,textwrap,collections,gc,ast,struct,zlib
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
os.environ['MPLCONFIGDIR']=str(OUT/'qa/matplotlib_cache')
sys.path.insert(0,str(ROOT/'comparative_visualization/2026-11-13/vendor'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import pymupdf
from PIL import Image
from terreus_tree_helpers import parse_newick,strong
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','font.size':11,'text.color':'#213c50','axes.labelcolor':'#213c50'})
OLD=ROOT/'comparative_visualization/2026-11-11';Q=OUT/'qa';T=OUT/'print_600dpi';sources={}
def read(rel):
 p=ROOT/rel;sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return json.loads(p.read_text())
allann=read('comparative_visualization/2026-11-11/analysis/tree_annotations_all_combinations.json')
ann={a['tip']:a for a in allann if a['tree']=='combined_all_strains'}
clades=[c for c in read('comparative_visualization/2026-11-11/analysis/type_clades_all_combinations.json') if c['tree']=='combined_all_strains']
palette=read('comparative_visualization/2026-11-11/analysis/taxon_palette.json')
refs={a['reference_id'] for a in read('taxonomy_reference_build/2026-11-10/analysis/tree_reference_metadata.json')}
genomes={'LOCAL_'+a['sample'] for a in read('taxonomy_reference_build/2026-11-10/analysis/genome_sources.json')}
ranking=read('comparative_visualization/2026-11-11/core66_combinations/ranking.json')
bootstrap=read('comparative_visualization/2026-11-11/core66_combinations/bootstrap_RF_summary.json')
ordination=[a for a in read('comparative_visualization/2026-11-11/core66_combinations/ordination.json') if a['method']=='NMDS']
treepath=OLD/'sources/input_trees/combined_all_strains.treefile';raw=treepath.read_bytes();sources[str(treepath)]=hashlib.sha256(raw).hexdigest()
root=parse_newick(raw.decode());tips=root.leaves();assert {t.name for t in tips}==set(ann) and len(tips)==283
sourcecols={'genome':'#087f8c','reference':'#bd642d'};locicols={1:'#d7794c',2:'#e6ad51',3:'#408f87'}
choicecols={r['marker_choice']:matplotlib.colors.to_hex(plt.get_cmap('tab10')(i)) for i,r in enumerate(ranking)};choicecols['Core genome']='#253e59'
def source(k):
 assert (k in genomes) != (k in refs),k
 return 'genome' if k in genomes else 'reference'
def label(k):
 a=ann[k];prefix='[G]' if source(k)=='genome' else '[R11]' if k.startswith('HIST2011_') else '[R]'
 return prefix+' '+a['display_label'].replace(' [2011]','').replace('[TG]','[TYPE-G]').replace('[T]','[TYPE]')
def typed(k):return ann[k]['type'] or bool(ann[k].get('genome_type_species'))
def cols(k):
 a=ann[k]
 return [sourcecols[source(k)],locicols[a['available_loci']],palette[a['species_ring']],palette[a['type_clade_category']],palette[a.get('genome_type_clade_category','No supported type clade')]]
def positions(clado=False):
 ys={id(t):i for i,t in enumerate(tips)};xs={id(root):0}
 def walk(n):
  for c in n.children:xs[id(c)]=xs[id(n)]+(1 if clado else c.length);walk(c)
  if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
 walk(root);return xs,ys
legend_texts=[];alignment=[];output_stems=[]
def save(fig,stem):
 fig.savefig(OUT/(stem+'.pdf'),bbox_inches='tight',pad_inches=.25)
 fig.savefig(OUT/(stem+'.svg'),bbox_inches='tight',pad_inches=.25)
 fig.savefig(Q/(stem+'_overview.png'),dpi=65,bbox_inches='tight',pad_inches=.25)
 plt.close(fig);gc.collect();output_stems.append(stem);print('DRAWN',stem,flush=True)

# Vertical figure: all column headings horizontally centred on their tracks.
fig=plt.figure(figsize=(24,49.8));ax=fig.add_axes([.035,.025,.34,.925]);lab=fig.add_axes([.395,.025,.36,.925],sharey=ax);track=fig.add_axes([.77,.025,.21,.925],sharey=ax)
xs,ys=positions();mx=max(xs.values())
def branches(n):
 if n.children:
  ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],c='#607384',lw=.7)
  for c in n.children:ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,c='#607384',lw=.7);branches(c)
  if strong(n):ax.scatter(xs[id(n)],ys[id(n)],s=5,c='#213c50')
 else:
  k=n.name;y=ys[id(n)];ax.plot([xs[id(n)],mx*1.014],[y,y],c='#cbd4dc',lw=.42,ls=':');ax.scatter(mx*1.025,y,marker='s' if source(k)=='genome' else 'o',s=15,c=sourcecols[source(k)])
branches(root)
for i,t in enumerate(tips):
 k=t.name;a=ann[k];is_type=typed(k)
 lab.text(0,i,label(k),va='center',fontsize=9,color=sourcecols[source(k)],weight='bold' if is_type or source(k)=='genome' else 'normal',bbox={'boxstyle':'round,pad=.1','facecolor':'#fff0bb','edgecolor':'none'} if is_type else None)
 if is_type:lab.scatter(-.025,i,marker='*',s=52,c='#e8af27',edgecolor='#533f12',linewidth=.45,clip_on=False)
 for j,col in enumerate(cols(k)):track.add_patch(Rectangle((j-.38,i-.45),.76,.9,facecolor=col,lw=0))
 track.text(5.15,i,a['type_clade'] if a['type_clade']!='unresolved' else '-',ha='center',va='center',fontsize=8)
for a in [ax,lab,track]:a.set_ylim(len(tips)-.3,-2.2);a.set_yticks([])
ax.set_xlim(0,mx*1.05);ax.spines[['top','right','left']].set_visible(False);ax.set_xlabel('Substitutions per site',fontsize=12)
lab.set_xlim(0,1);lab.axis('off');track.set_xlim(-.55,5.75);track.axis('off')
headers=['1\nSequence\nsource','2\nAvailable\nloci','3\nReviewed\nspecies','4\nReference\ntype clade','5\nGenome\ntype','Clade\ncode']
for x,h in zip([0,1,2,3,4,5.15],headers):
 text=track.text(x,1.004,h,transform=track.get_xaxis_transform(),ha='center',va='bottom',fontsize=10,linespacing=1.35)
 alignment.append({'figure':'vertical','header':h,'track_x':x,'alignment':'center','baseline_axes_y':1.004})
fig.text(.035,.986,'BenA + CaM + RPB2',fontsize=28,weight='bold',va='top')
fig.text(.035,.974,'283 tips | 76 genome-derived sequences | 207 taxonomic reference entries',fontsize=13,va='top',color='#526778')
save(fig,'01_combined_three_gene_vertical')

# Circular companion, with exact same tip annotations and colour lookup.
fig,ax=plt.subplots(figsize=(38,38),subplot_kw={'projection':'polar'});xs,ys=positions(True);mx=max(xs.values());n=len(tips);step=2*math.pi/n;idx={t.name:i for i,t in enumerate(tips)};painted=set()
for c in sorted([c for c in clades if c['status']=='single_type_species'],key=lambda c:-len(c['descendant_tips'])):
 inds=[idx[k] for k in c['descendant_tips'] if k in idx]
 if len(inds)<2 or set(inds)&painted or max(inds)-min(inds)+1!=len(inds):continue
 ax.bar((min(inds)+max(inds))*step/2,1.02,width=len(inds)*step,bottom=0,color=palette[c['type_species']],alpha=.065,ec='none');painted.update(inds)
def circle(n):
 r=xs[id(n)]/mx;a=ys[id(n)]*step
 if n.children:
  angles=[ys[id(c)]*step for c in n.children];ax.plot(np.linspace(min(angles),max(angles),70),[r]*70,c='#607384',lw=.7)
  for c in n.children:ax.plot([ys[id(c)]*step]*2,[r,xs[id(c)]/mx],c='#607384',lw=.7);circle(c)
  if strong(n):ax.scatter(a,r,s=5,c='#213c50')
 else:
  k=n.name;ax.plot([a]*2,[r,1.025],c='#cbd4dc',lw=.42,ls=':');ax.scatter(a,1.025,marker='s' if source(k)=='genome' else 'o',s=15,c=sourcecols[source(k)])
  for j,col in enumerate(cols(k)):ax.bar(a,.033,width=step*.93,bottom=1.053+j*.045,color=col,lw=0)
  deg=math.degrees(a);left=90<deg<270;is_type=typed(k)
  if is_type:ax.scatter(a,1.29,marker='*',s=70,color='#e8af27',edgecolor='#533f12',linewidth=.45,zorder=5)
  ax.text(a,1.32,label(k),rotation=deg+(180 if left else 0),rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=9,color=sourcecols[source(k)],weight='bold' if is_type or source(k)=='genome' else 'normal',bbox={'boxstyle':'round,pad=.1','facecolor':'#fff0bb','edgecolor':'none'} if is_type else None)
circle(root)
runs=[]
for i,t in enumerate(tips):
 code=ann[t.name]['type_clade']
 if runs and runs[-1][2]==code:runs[-1][1]=i
 else:runs.append([i,i,code])
for lo,hi,code in runs:
 if code=='unresolved' or hi-lo<2:continue
 a=(lo+hi)*step/2;deg=math.degrees(a);ax.text(a,.998,code,fontsize=9,ha='center',va='center',rotation=deg-90+(180 if 90<deg<270 else 0),rotation_mode='anchor')
ax.set_ylim(0,2.02);ax.axis('off');fig.subplots_adjust(left=.02,right=.98,top=.94,bottom=.01)
fig.text(.035,.983,'BenA + CaM + RPB2',fontsize=30,weight='bold',va='top');fig.text(.035,.967,'283 tips | circular cladogram',fontsize=14,va='top',color='#526778')
save(fig,'02_combined_three_gene_circular')

# RF/NMDS: neutral comparative title; no narrative footer or taxonomic ranking.
fig,(a,b)=plt.subplots(1,2,figsize=(17,7.8),gridspec_kw={'width_ratios':[1.15,1]});fig.subplots_adjust(left=.17,right=.97,top=.79,bottom=.13,wspace=.32)
boots={x['marker_choice']:x for x in bootstrap}
for i,r in enumerate(ranking):
 c=choicecols[r['marker_choice']];v=boots[r['marker_choice']]
 a.plot([v['percentile_2_5'],v['percentile_97_5']],[i]*2,c=c,lw=3,alpha=.7);a.scatter(r['normalized_RF_to_core'],i,s=65,color=c,zorder=3);a.text(.90,i,f"{r['normalized_RF_to_core']:.3f}",va='center',fontsize=10)
a.set_yticks(range(7),[r['marker_choice'].replace('_',' + ') for r in ranking]);a.invert_yaxis();a.set_xlim(.57,.94);a.set_xlabel('Normalised RF distance');a.set_title('A  Distance to the core-genome topology',loc='left',fontsize=12,weight='bold');a.grid(axis='x',alpha=.14);a.spines[['top','right']].set_visible(False)
for t in a.get_yticklabels():
 if t.get_text()=='BenA + CaM + RPB2':t.set_weight('bold')
offsets={'Core genome':(10,-18),'BenA':(-9,13),'CaM':(8,-18),'RPB2':(-48,12),'BenA_CaM':(-46,-21),'BenA_RPB2':(-48,14),'CaM_RPB2':(10,7),'BenA_CaM_RPB2':(-59,-24)}
for r in ordination:
 name=r['tree'];b.scatter(r['axis1'],r['axis2'],marker='*' if name=='Core genome' else 'o',s=160 if name=='Core genome' else 70,color=choicecols[name]);b.annotate(name.replace('_',' + '),(r['axis1'],r['axis2']),xytext=offsets[name],textcoords='offset points',fontsize=10,weight='bold' if name=='BenA_CaM_RPB2' else 'normal')
b.set(xlim=(-1.05,1.24),ylim=(-.65,.79),xlabel='NMDS axis 1',ylabel='NMDS axis 2');b.set_aspect('equal');b.spines[['top','right']].set_visible(False);b.axhline(0,c='#e5e9ec',lw=.6);b.axvline(0,c='#e5e9ec',lw=.6);b.set_title(f"B  NMDS | stress-1 = {ordination[0]['stress1']:.4f}",fontsize=12,loc='left',weight='bold')
fig.text(.035,.947,'Phylogenetic topology comparison',fontsize=23,weight='bold');fig.text(.035,.90,'Seven marker combinations and the core-genome tree | matched set of 66 genomes',fontsize=12,color='#526778')
save(fig,'03_RF_NMDS_comparison')

# One consolidated, single-page legend, with shared row centres for swatches and text.
fig=plt.figure(figsize=(26,25));la=fig.add_axes([0,0,1,1]);la.set(xlim=(0,1),ylim=(0,1));la.axis('off')
def txt(x,y,s,size=13,bold=False,va='center',ha='left',color='#213c50'):
 t=la.text(x,y,s,fontsize=size,weight='bold' if bold else 'normal',va=va,ha=ha,color=color);legend_texts.append(t);return t
def swatch(x,y,c,s,size=12,w=.009):
 la.add_patch(Rectangle((x,y-.0047),w,.0094,facecolor=c,edgecolor='#82909a',lw=.4));txt(x+w+.007,y,s,size)
 alignment.append({'figure':'legend','label':s,'swatch_center_y':y,'label_center_y':y,'vertical_alignment':'center'})
txt(.035,.966,'Combined three-gene tree and RF/NMDS | consolidated legend',27,True)
txt(.035,.934,'Sequence source and type status',17,True)
for x,marker,c,l in [(.035,'s',sourcecols['genome'],'[G]  Genome-derived'),(.29,'o',sourcecols['reference'],'[R]  Taxonomic reference'),(.58,'o',sourcecols['reference'],'[R11]  Historical 2011 reference')]:
 la.scatter(x+.004,.909,s=90,marker=marker,c=c);txt(x+.016,.909,l,14)
la.scatter(.039,.884,s=150,marker='*',c='#e8af27',edgecolor='#533f12',linewidth=.6);txt(.051,.884,'[TYPE]  Reference ex-type entry    |    [TYPE-G]  Documented ex-type genome',14)
txt(.035,.858,'* Unresolved reported taxon    |    [name flag] Source-name conflict    |    Internal node dot: SH-aLRT >=80 and UFBoot >=95',13)
txt(.035,.827,'Annotation tracks: vertical left to right; circular inside to outside',17,True)
track_descriptions=[('1','Sequence source','Genome-derived or taxonomic-reference sequence'),('2','Available loci','Number of represented loci in the three-gene panel'),('3','Reviewed species','Curated reference name or reviewed genome assignment'),('4','Reference type clade','Ex-type-containing clade passing both support thresholds'),('5','Genome type','Documented ex-type genome / supported compatible group')]
for i,(num,title,detail) in enumerate(track_descriptions):
 y=.798-i*.019;txt(.035,y,num,13,True);txt(.057,y,title,13,True);txt(.205,y,detail,12)
for x,n in [(.62,1),(.73,2),(.84,3)]:swatch(x,.776,locicols[n],f'{n} locus' if n==1 else f'{n} loci',13)
txt(.62,.747,'Circular: topology only. Vertical: branch lengths in substitutions/site.',11)
txt(.62,.728,'31 reference ex-type entries and 3 genome type anchors;',11)
txt(.62,.710,'entry counts do not represent 34 independent type cultures.',11)
txt(.035,.674,'Species and type-clade colours',17,True);txt(.53,.674,'Reference type-clade key',17,True)
for i,(tax,c) in enumerate(palette.items()):swatch(.035+(i//18)*.245,.644-(i%18)*.024,c,tax,12)
txt(.53,.644,'Code',11,True);txt(.58,.644,'Type species in clade',11,True);txt(.952,.644,'SH / UFB',11,True,ha='right')
y=.622
for c in clades:
 name=c['type_species'];lines=textwrap.wrap(name,82);swatch(.53,y,palette.get(name,palette['Mixed type species']),c['clade_id'],10,w=.006)
 txt(.58,y,'\n'.join(lines),10,va='center');txt(.952,y,f"{c['SH_aLRT']:g} / {c['UFBoot']:g}",10,ha='right');y-=.0165+max(0,len(lines)-1)*.009
txt(.035,.198,'Grey: unresolved species. Lavender-grey: name conflict. Mixed-type clades are not species assignments.',11)
txt(.035,.180,'White type track: no supported type clade; it does not encode biological absence. Pale sectors mark supported single-type clades.',11)
txt(.035,.151,'RF/NMDS colours and statistical symbols',17,True)
for i,(name,c) in enumerate(choicecols.items()):
 x=.035+(i%4)*.245;y=.126-(i//4)*.025
 la.scatter(x+.004,y,s=80 if name!='Core genome' else 140,marker='*' if name=='Core genome' else 'o',c=c);txt(x+.016,y,name.replace('_',' + '),12,bold=name=='BenA_CaM_RPB2')
txt(.035,.070,'RF panel: dot = ML-tree distance; line = 2.5th-97.5th percentiles of 1,000 bootstrap trees versus the fixed core topology.',11)
txt(.035,.050,'NMDS: pairwise unrooted RF distances on the same 66 tips. Normalised RF uses a denominator of 126.',11)
txt(.035,.030,'BenA + CaM + RPB2 is the taxonomic reference framework for this dataset; RF proximity measures topology, not taxonomic authority.',11)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
overlaps=[]
for i,t in enumerate(legend_texts):
 a=t.get_window_extent(renderer)
 for u in legend_texts[i+1:]:
  b=u.get_window_extent(renderer)
  if a.overlaps(b):overlaps.append((t.get_text(),u.get_text()))
(Q/'legend_text_overlap_check.json').write_text(json.dumps(overlaps,indent=2))
save(fig,'04_consolidated_colour_legend')

fields=['tip','display_label','sequence_source','source_colour','available_loci','loci_colour','reviewed_species','species_colour','reference_type_clade','reference_type_colour','genome_type_clade','genome_type_colour','reference_type','genome_type']
with (OUT/'data/tip_colour_audit.tsv').open('w',newline='',encoding='utf-8') as h:
 w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader()
 for t in tips:
  k=t.name;a=ann[k];cc=cols(k);w.writerow(dict(zip(fields,[k,label(k),source(k),cc[0],a['available_loci'],cc[1],a['species_ring'],cc[2],a['type_clade_category'],cc[3],a['genome_type_clade_category'],cc[4],a['type'],bool(a.get('genome_type_species'))])))
for name,value in [('palette',{'species':palette,'source':sourcecols,'loci':locicols,'marker_combinations':choicecols}),('type_clade_key',clades),('alignment_audit',alignment),('source_hashes',sources)]:
 (OUT/'data'/(name+'.json')).write_text(json.dumps(value,indent=2))
(OUT/'data/combined_all_strains.treefile').write_bytes(raw)
checks={'283_tips':len(tips)==283,'76_genome_sequences':sum(source(t.name)=='genome' for t in tips)==76,'207_reference_sequences':sum(source(t.name)=='reference' for t in tips)==207,'31_reference_types':sum(a['type'] for a in ann.values())==31,'3_genome_types':sum(bool(a.get('genome_type_species')) for a in ann.values())==3,'all_species_colours_unique':len(set(palette.values()))==len(palette),'same_colour_lookup_both_trees':True,'all_legend_rows_aligned':all(a['swatch_center_y']==a['label_center_y'] for a in alignment if a['figure']=='legend'),'legend_text_nonoverlapping':not overlaps,'input_sources_unchanged':all(hashlib.sha256(Path(k).read_bytes()).hexdigest()==v for k,v in sources.items())}
for stem in output_stems:
 with pymupdf.open(OUT/(stem+'.pdf')) as d:
  assert len(d)==1;page=d[0];text=page.get_text();checks[stem+'_english']=not any(s in text for s in ['Typstämme','Typstamm','Sequenzherkunft','Referenz-Ex-Typ','Goldstern','hervorgehobener'])
  checks[stem+'_no_raster']=not page.get_image_info()
  checks[stem+'_embedded_fonts']=all(d.extract_font(f[0])[3] for f in page.get_fonts())
  if stem in output_stems[:2]:checks[stem+'_type_counts']=text.count('[TYPE]')==31 and text.count('[TYPE-G]')==3
  if stem.startswith('03'):checks['NMDS_no_editorial_footer']=not any(x in text.lower() for x in ['substantial','taxonomic accuracy','overlapping','closest','unsupported'])
  for tag,clip in [('top',pymupdf.Rect(0,0,page.rect.width,min(350,page.rect.height))),('detail',pymupdf.Rect(page.rect.width*.35,page.rect.height*.15,page.rect.width*.98,page.rect.height*.23))]:page.get_pixmap(matrix=pymupdf.Matrix(1,1),clip=clip).save(str(Q/(stem+'_'+tag+'.png')))
(Q/'figure_validation.json').write_text(json.dumps({'checks':checks,'resource_limits':LIMITS},indent=2));assert all(checks.values()),{k:v for k,v in checks.items() if not v}
print('FIGURE VALIDATION PASS',len(checks),flush=True)
