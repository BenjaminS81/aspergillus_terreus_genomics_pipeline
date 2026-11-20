"""Single-page vector figures, 300-dpi images and a linked scientific review."""
from pathlib import Path
import os,json,sys,gc,html,shutil,hashlib,textwrap
ROOT=Path(__file__).resolve().parent
from circos_review_common import constrain_process
LIMITS=constrain_process()
from terreus_tree_helpers import parse_newick
OUT=ROOT/'stakeholder_release/2026-11-20_v4_clade_congruence_cgmlst'
os.environ['MPLCONFIGDIR']=str(OUT/'qa/mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
import fitz
def read(n):return json.loads((OUT/'data'/n).read_text())
MAN=read('analysis_manifest.json');KEYS=MAN['tree_keys'][1:];ORDER=MAN['core_order']
ANN={a['tip']:a for a in read('core_tip_assignments.json')};CL=read('clade_congruence.json');G=read('global_comparisons.json');INV=read('cgmlst_candidate_inventory.json')
PAL=json.loads((ROOT/'stakeholder_release/2026-11-17_v3_taxonomy_review/data/palette.json').read_text())['species']
TITLE={k:k.replace('_',' + ') for k in KEYS}
TITLE['BenA_CaM_RPB2']='BenA + CaM + RPB2 | inferred on 66 genomes'
TITLE['Three_gene_reference230']='BenA + CaM + RPB2 | reference tree reduced from 230 to 66'
SHORT={k:k.replace('_',' + ') for k in KEYS};SHORT[KEYS[-2]]='3 genes\n66 inferred';SHORT[KEYS[-1]]='3 genes\n230 reduced'
INK='#203849';GRAY='#70818f';GOLD='#d9a72f';YES='#c9e7e4';NO='#f5d4b9';NA='#edf0f3'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,'pdf.fonttype':42,'svg.fonttype':'none'})
OUTPUTS=[];QA=[]
def output(fig,stem):
 if '--legend-only' in sys.argv and stem!='12_shared_legend':
  plt.close(fig);OUTPUTS.append(stem)
  QA.append(next(a for a in json.loads((OUT/'qa/figure_checks.json').read_text()) if a['figure']==stem));return
 p=OUT/'figures'/stem
 fig.savefig(p.with_suffix('.pdf'),facecolor='white');fig.savefig(p.with_suffix('.svg'),facecolor='white');plt.close(fig);gc.collect()
 with fitz.open(p.with_suffix('.pdf')) as doc:
  assert len(doc)==1
  page=doc[0];pix=page.get_pixmap(dpi=300);pix.save(str(p.with_suffix('.png')));size=[pix.width,pix.height];del pix
  page.get_pixmap(dpi=55).save(str(OUT/'qa'/(stem+'.png')))
  oob=[]
  for b in page.get_text('dict')['blocks']:
   for l in b.get('lines',[]):
    for s in l['spans']:
     x0,y0,x1,y1=s['bbox']
     if x0<-.1 or y0<-.1 or x1>page.rect.width+.1 or y1>page.rect.height+.1:oob.append(s['text'])
  QA.append({'figure':stem,'pages':1,'raster_dpi':300,'raster_pixels':size,'out_of_page_text':oob})
 OUTPUTS.append(stem);print('FIGURE',stem,flush=True);gc.collect()
def positions(r):
 tips=r.leaves();ys={id(t):i for i,t in enumerate(tips)};xs={id(r):0}
 def walk(n):
  for c in n.children:xs[id(c)]=xs[id(n)]+c.length;walk(c)
  if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
 walk(r);return tips,xs,ys
TREES={k:parse_newick((OUT/'data'/(k+'.matched66.treefile')).read_text()) for k in MAN['tree_keys']}
def branchplot(ax,r,xs,ys,extent):
 def draw(n):
  if n.children:
   ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],color=GRAY,lw=.75)
   for c in n.children:
    species={ANN[t.name]['assigned_species'] for t in c.leaves()};co=PAL[next(iter(species))] if len(species)==1 else GRAY
    ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,c=co,lw=1.1);draw(c)
  else:ax.plot([xs[id(n)],extent],[ys[id(n)]]*2,c='#cbd2d8',lw=.45,ls=':')
 draw(r);ax.set_xlim(0,extent);ax.set_ylim(65.65,-.65);ax.set_yticks([])
 for s in ('left','right','top'):ax.spines[s].set_visible(False)
 ax.spines['bottom'].set_color('#b9c3ca');ax.tick_params(axis='x',labelsize=9);ax.set_xlabel('Substitutions per site',fontsize=10)
for num,k in enumerate(KEYS,1):
 core=TREES['Core'];marker=TREES[k];lt,lx,ly=positions(core);rt,rx,ry=positions(marker)
 extent=max(max(lx.values()),max(rx.values()))*1.04
 fig=plt.figure(figsize=(30,19));left=fig.add_axes([.025,.065,.155,.86]);right=fig.add_axes([.82,.065,.155,.86])
 labL=fig.add_axes([.187,.065,.246,.86]);links=fig.add_axes([.435,.065,.13,.86]);labR=fig.add_axes([.57,.065,.245,.86])
 branchplot(left,core,lx,ly,extent);branchplot(right,marker,rx,ry,extent);right.invert_xaxis()
 for ax in (labL,labR,links):ax.set_ylim(65.65,-.65);ax.set_xlim(0,1);ax.axis('off')
 yL={t.name:ly[id(t)] for t in lt};yR={t.name:ry[id(t)] for t in rt}
 for t in ORDER:
  a=ANN[t];co=PAL[a['assigned_species']];links.plot([0,1],[yL[t],yR[t]],color=co,alpha=.65,lw=.85)
  for ax,y in ((labL,yL[t]),(labR,yR[t])):
   ax.add_patch(Rectangle((.001,y-.34),.009,.68,color=co,lw=0))
   ax.text(.032,y,a['display_label'],fontsize=8.4,va='center',weight='bold' if a['is_type'] else 'normal')
   if a['is_type']:ax.scatter(.02,y,marker='*',s=42,c=GOLD,edgecolors='#71551b',lw=.4)
 fig.text(.025,.974,'Core genome vs '+TITLE[k],fontsize=23,weight='bold')
 fig.text(.025,.947,'66 identical genome labels on both sides | fixed core order | branch-length scale shared within this pair',fontsize=12,color=GRAY)
 left.set_title('CORE GENOME',fontsize=13,loc='left',pad=17);right.set_title('MARKER TREE',fontsize=13,loc='right',pad=17)
 output(fig,f'{num:02d}_tanglegram_core_vs_{k}')
# Clade boundary recovery: color describes the point estimate; text is projected bootstrap frequency.
groups=list(dict.fromkeys(a['group'] for a in CL));trees=MAN['tree_keys'];lookup={(a['group'],a['tree']):a for a in CL}
fig,ax=plt.subplots(figsize=(17,9));fig.subplots_adjust(left=.20,right=.98,bottom=.21,top=.82)
vals=np.array([[-1 if lookup[g,k]['strict_split_recovered'] is None else int(lookup[g,k]['strict_split_recovered']) for k in trees] for g in groups])
ax.imshow(vals,cmap=ListedColormap([NA,NO,YES]),vmin=-1,vmax=1,aspect='auto')
for i,g in enumerate(groups):
 for j,k in enumerate(trees):
  a=lookup[g,k];f=a['bootstrap_split_frequency_pct'];txt='N/A' if vals[i,j]<0 else ('Yes' if vals[i,j] else 'No') if k=='Core' else f'{f:.1f}%'
  ax.text(j,i,txt,ha='center',va='center',fontsize=11)
ax.set_xticks(range(len(trees)),['Core']+[SHORT[k] for k in KEYS],rotation=30,ha='right');ax.set_yticks(range(len(groups)),[g+f"  (n={lookup[g,'Core']['n']})" for g in groups]);ax.tick_params(length=0)
fig.text(.035,.945,'Recovery of type-derived genome groups',fontsize=23,weight='bold')
fig.text(.035,.892,'Cell color: exact group recovered / not recovered | values: frequency among 1,000 projected UFBoot trees',fontsize=11)
for x,col,txt in ((.21,YES,'Recovered'),(.40,NO,'Not recovered'),(.62,NA,'Singleton: not assessable')):
 fig.add_artist(Rectangle((x,.045),.013,.019,transform=fig.transFigure,color=col));fig.text(x+.02,.05,txt,fontsize=11)
output(fig,'09_clade_boundary_recovery')
validgroups=[g for g in groups if lookup[g,'Core']['n']>=4]
mat=np.array([[lookup[g,k]['within_group_normalized_RF_to_core'] for k in KEYS] for g in validgroups])
fig,ax=plt.subplots(figsize=(16,8));fig.subplots_adjust(left=.22,right=.89,bottom=.24,top=.81)
im=ax.imshow(mat,vmin=0,vmax=1,cmap='magma_r',aspect='auto')
for i in range(len(validgroups)):
 for j in range(len(KEYS)):ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',color='white' if mat[i,j]>.52 else INK,fontsize=13)
ax.set_xticks(range(len(KEYS)),[SHORT[k] for k in KEYS],rotation=30,ha='right');ax.set_yticks(range(len(validgroups)),[g+f"  (n={lookup[g,'Core']['n']})" for g in validgroups]);ax.tick_params(length=0)
cax=fig.add_axes([.92,.24,.016,.57]);fig.colorbar(im,cax=cax,label='Normalized RF')
fig.text(.035,.94,'Relationships within the same species-assigned set',fontsize=22,weight='bold')
fig.text(.035,.88,'Each cell compares the same members in the core and marker tree; species-boundary recovery is assessed separately.',fontsize=11)
output(fig,'10_within_group_topology_RF')
tips=read('tip_quartet_discordance.json');tmap={(a['tip'],a['tree']):a for a in tips};mat=np.array([[tmap[t,k]['discordance_pct'] for k in KEYS] for t in ORDER])
fig,ax=plt.subplots(figsize=(17,21));fig.subplots_adjust(left=.37,right=.89,bottom=.085,top=.925)
im=ax.imshow(mat,vmin=0,vmax=100,cmap='magma_r',aspect='auto')
ax.set_yticks(range(66),[ANN[t]['display_label'] for t in ORDER],fontsize=7.3);ax.set_xticks(range(8),[SHORT[k] for k in KEYS],rotation=25,ha='right',fontsize=10);ax.tick_params(length=0)
for i in range(66):
 for j in range(8):ax.text(j,i,f'{mat[i,j]:.0f}',ha='center',va='center',fontsize=8,color='white' if mat[i,j]>52 else INK)
cax=fig.add_axes([.925,.28,.015,.43]);fig.colorbar(im,cax=cax,label='Discordant resolved quartets containing this genome (%)')
fig.text(.035,.973,'Where does the topological disagreement occur?',fontsize=23,weight='bold')
fig.text(.035,.947,'All 720,720 quartets evaluated per comparison | same 66 genomes | core order retained',fontsize=11)
output(fig,'11_genome_quartet_discordance')
# One shared legend. Keep interpretation rules here, not prose on each tanglegram.
fig=plt.figure(figsize=(16,11));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off');fig.text(.045,.94,'Shared figure legend',fontsize=26,weight='bold')
fig.text(.045,.886,'Colors identify inherited three-gene assignments',fontsize=15,weight='bold')
species=sorted({a['assigned_species'] for a in ANN.values()})
for i,sp in enumerate(species):
 y=.836-i*.056;ax.add_patch(Rectangle((.05,y-.010),.019,.025,color=PAL[sp]));members=[a for a in ANN.values() if a['assigned_species']==sp]
 culture=next((a['display_label'].split('TYPE: ')[1] for a in members if a['is_type']),'; '.join(members[0]['type_cultures']))
 fig.text(.082,y,sp,fontsize=13);fig.text(.29,y,f"n={len(members)}",fontsize=11);fig.text(.345,y,culture,fontsize=11)
fig.text(.045,.286,'Actual type genomes shown in the matched trees',fontsize=14,weight='bold')
for i,a in enumerate([a for a in ANN.values() if a['is_type']]):
 y=.246-i*.048;ax.scatter(.058,y+.004,marker='*',s=100,c=GOLD,edgecolor='#71551b',lw=.5);fig.text(.082,y,a['sample']+' | '+a['assigned_species']+' | '+a['display_label'].split('TYPE: ')[1],fontsize=10)
fig.text(.55,.886,'How to read the comparisons',fontsize=15,weight='bold')
paras=[
 ('Leaf labels','Inferred species (original source taxon) | original strain ID. TYPE labels and gold stars identify only the three documented type genomes actually present.'),
 ('Connecting lines','Each line joins the same genome. Species colors are fixed across trees. Marker colors do not constitute independent taxonomic assignments.'),
 ('Linear scale','Branch lengths are substitutions per site; each pair shares a numerical scale. Trees use different fitted models, so lengths are not elapsed time. The floccosus tip is a display root only.'),
 ('Clade boundary heatmap','Green: the exact member set is separated by a split. Orange: it is not. Values show its frequency in 1,000 restricted UFBoot trees. Gray: singleton, not a multi-genome test.'),
 ('Within-group RF','0 = identical unrooted splits; 1 = no shared internal splits. The member set is fixed even when it is not monophyletic. Short or weak branches can increase RF.'),
 ('Genome quartet heatmap','Percentage of resolved quartets containing that genome whose topology differs from the core. This locates disagreement; it is not a probability that a genome is misidentified.'),
 ('Support and taxonomic scope','Core bootstrap support is unavailable. Restricted bootstrap frequencies are not full type-clade support. Groups defined by the three-gene tree are not independently validated by recovery in that same tree.')]
y=.83
for title,txt in paras:
 fig.text(.55,y,title,fontsize=12,weight='bold');wrapped=textwrap.fill(txt,76);fig.text(.55,y-.026,wrapped,fontsize=10,va='top',linespacing=1.4);y-=.026+.0195*(wrapped.count('\n')+1)+.025
output(fig,'12_shared_legend')
# Scientific report and explicit future gates, with primary documentation links.
GATES=[
 {'gate':'C0','status':'PARTIAL','input':'Assembly accessions, sample/culture keys, current BUSCO and assembly QC; type references','action':'Resolve duplicate sequence hashes versus biological replicates; include missing public genomes only after QC. Obtain verified type genomes for terreus and montoensis; broaden section-wide sampling.','output':'Frozen training/holdout panel and provenance manifest; explicit exclusions.','criterion':'Unique sequence identity and culture identity separately tracked; no confident species scheme for unresolved assignments.'},
 {'gate':'C1','status':'OPEN','input':'Genome FASTA + harmonized fungal GFF3/CDS/protein files','action':'Validate coordinate-consistent exon-joined CDS, translation, genetic code, internal stops and one representative isoform per gene. Distinguish missing assembly sequence from failed annotation.','output':'Versioned fungal CDS and protein catalog with per-locus callability codes.','criterion':'Fungal annotation pipeline validated across assembly sources; protein-only/transcript-only data cannot by themselves supply comparable nucleotide genome alleles.'},
 {'gate':'C2','status':'SCREENED, NOT VALIDATED','input':'14,964 current orthogroups; 4,661 strict single-copy candidates','action':'Review orthology with protein trees, duplication histories, reciprocal matches and synteny. Screen paralogs, repeats and unstable models. Re-evaluate on a species-balanced independent panel.','output':'Candidate loci with stable IDs, coordinates, copy-state evidence and exclusions.','criterion':'Exactly one defensible ortholog per callable genome; copy count alone is insufficient.'},
 {'gate':'C3','status':'OPEN','input':'Validated loci and representative diversity panel','action':'Pilot exon-joined nucleotide CDS alleles using chewBBACA CDS mode or a validated allele catalog/caller. Benchmark intron handling against curated loci. Keep intronic genomic haplotypes as a separately named optional analysis.','output':'Reproducible allele sequences, hashes and immutable allele IDs; schema version and allele profiles.','criterion':'No arbitrary paralog selection; missing, ambiguous, truncated and multicopy calls remain distinct from valid allele IDs.'},
 {'gate':'C4','status':'DESIGN READY','input':'Validated section-wide catalog and defensible species groups','action':'Define a stable section-wide backbone plus species-specific extensions. Retain global locus identifiers and version membership lists. Compare 95%, 99% and 100% occupancy thresholds on training data.','output':'Section backbone and nested species views/extensions; no silently changing core.','criterion':'Thresholds are proposed pilot settings, selected using holdout callability and stability, not declared universal fungal standards.'},
 {'gate':'C5','status':'OPEN','input':'Independent holdout genomes, alternate assemblies/platforms and technical replicates','action':'Measure allele-call reproducibility, dropout, discrimination, saturation and agreement with validated core SNP/phylogenomic relationships. Test annotation-version sensitivity and species imbalance.','output':'Validation report, failure matrix and minimum-callability policy.','criterion':'Report differences/shared callable loci per pair; no missing-as-zero distance, no universal outbreak cutoff without empirical validation.'},
 {'gate':'C6','status':'OPEN','input':'Frozen allele catalog and validation evidence','action':'Package allele nomenclature, schema changelog, test panel and reproducible caller. Keep lovastatin architecture/phenotypes as accessory traits for structure-aware evolutionary analyses.','output':'Research cgMLST release and reproducible profiles, with optional later public registry submission.','criterion':'Research pilot clearly separated from externally validated nomenclature; ecological or clinical thresholds require independent evidence.'}
]
(OUT/'data/cgmlst_gates.json').write_text(json.dumps(GATES,indent=2))
def esc(x):return html.escape(str(x))
def table(rows,columns):
 return '<table><thead><tr>'+''.join('<th>'+esc(c[0])+'</th>' for c in columns)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(c[1](r))+'</td>' for c in columns)+'</tr>' for r in rows)+'</tbody></table>'
gallery=''.join(f'<article><h3>{esc(s.replace("_"," "))}</h3><a href="figures/{s}.pdf">Single-page vector PDF</a> | <a href="figures/{s}.png">300-dpi PNG</a> | <a href="figures/{s}.svg">SVG</a><a href="figures/{s}.pdf"><img loading="lazy" src="qa/{s}.png"></a></article>' for s in OUTPUTS)
htmltext='''<!doctype html><html><head><meta charset="utf-8"><title>Terrei clade congruence and cgMLST assessment</title><style>body{font:16px/1.55 system-ui,sans-serif;color:#203849;max-width:1400px;margin:36px auto;padding:0 24px}h1,h2{color:#173e53}table{border-collapse:collapse;width:100%;font-size:13px;margin:20px 0}th,td{padding:9px;border:1px solid #d9e1e7;text-align:left;vertical-align:top}th{background:#eaf1f5}article{border:1px solid #d9e1e7;padding:18px;margin:25px 0}img{display:block;max-width:100%;max-height:1100px;margin:12px auto}.note{padding:18px;background:#f4f7f9}a{color:#126a8b}code{white-space:pre-wrap}</style></head><body>'''
htmltext+='''<h1>Section Terrei: clade-level congruence and cgMLST feasibility</h1><p>Version 2026-11-20_v4. New analysis of existing trees and copy-count data; source files and running inference/QC jobs were not modified.</p><div class="note"><b>Scope:</b> all eight comparisons contain exactly the same 66 genome labels. Seven marker combinations were inferred on those 66 genomes. The eighth is the current 230-tip taxonomic three-gene tree restricted to the same cohort. All 66 aligned sequence rows match exactly across these input versions (2,076 columns); taxon sampling and inference differ. Figures are separate one-page PDFs, with vector SVG and 300-dpi PNG companions.</div>'''
htmltext+='<h2>What was added</h2><ol><li>Exact group recovery and projected bootstrap frequencies per species-assigned group.</li><li>Within-group RF, keeping the member set identical even when the group is not a core clade.</li><li>Exhaustive comparison of all 720,720 four-genome subsets, with disagreement localized per genome.</li><li>Core-edge recovery and marker conflict audit; nearest-split precision/recall and explicit extra/missing members.</li><li>Eight branch-length-preserving linear tanglegrams with fixed core order and deterministic marker rotations.</li><li>A copy-count feasibility screen for a section-wide cgMLST backbone and species extensions.</li></ol>'
htmltext+='<h2>Global results</h2>'+table(G,[('Marker comparison',lambda r:TITLE[r['tree']]),('Normalized RF',lambda r:f"{r['normalized_RF']:.3f}"),('Shared core splits / 63',lambda r:r['shared_core_splits']),('Quartet agreement',lambda r:f"{r['quartet_agreement_pct']:.2f}%"),('RF, edges >1.1e-6',lambda r:f"{r['normalized_RF_edges_over_1_1e_6']:.3f}")])
htmltext+='''<p>The de novo three-gene tree has the lowest RF, whereas CaM has the highest quartet agreement. These metrics weight disagreements differently, so neither ranking establishes a universally best marker or taxonomic truth. The reference-rich three-gene tree remains the type-context analysis; its reduced tree is not equivalent to re-inference on 66 genomes.</p>'''
focus=[a for a in CL if a['tree']==KEYS[-1] and a['n']>=2]
htmltext+='<h2>Current reference-rich three-gene groups against the core</h2>'+table(focus,[('Group',lambda a:a['group']),('n',lambda a:a['n']),('Core split recovered',lambda a:lookup[a['group'],'Core']['strict_split_recovered']),('Three-gene restricted split frequency',lambda a:f"{a['bootstrap_split_frequency_pct']:.1f}%"),('Within-set normalized RF',lambda a:'N/A (<4 genomes)' if a['within_group_normalized_RF_to_core'] is None else f"{a['within_group_normalized_RF_to_core']:.3f}")])
htmltext+='''<p>Terreus and montoensis are individually non-monophyletic in the current core topology; their 40-genome union is recovered. Neither corresponding type genome is present in this mixed core region. Citrinoterreus, hortae and pseudoterreus retain their boundaries but disagree in internal strain relationships. These are conditional comparisons to a core topology without comparable bootstrap evidence. In particular, a high restricted genome-group bootstrap is not the support of the enclosing full reference/type clade: pruning reference tips can remove the competing placements responsible for low full-tree support.</p><p>Species groups were defined from the reference-rich three-gene tree. Their recovery in that tree is therefore expected by construction, not independent validation. Singletons cannot test monophyly; two-genome groups can test a boundary but cannot support a within-group RF analysis. Several labels refer to identical assembly sequences, so the 66-tip cohort is not necessarily 66 independent biological replicates.</p>'''
htmltext+='<h2>cgMLST: a credible starting point, not yet a validated scheme</h2><p>The existing matrix contains 14,964 orthogroups. Exactly 4,661 are single-copy across all 66 labels; this is precisely the existing core-locus list. Allowing up to 5% absence while forbidding duplicate models yields 6,129 candidates. The current panel contains nine assigned species and 64 distinct assembly sequence hashes; it is not section-wide coverage. These are gene-model counts, not validated, portable nucleotide allele definitions.</p>'
htmltext+=table(INV['panels'],[('Panel',lambda r:r['panel']),('Genome labels',lambda r:r['tip_labels']),('Distinct sequence hashes',lambda r:r['distinct_sequence_hashes_or_unresolved_ids']),('100% single-copy candidates',lambda r:r['all_single_copy_100pct']),('95% occupancy; no duplicates',lambda r:r['at_least_95pct_single_copy_no_duplicates'])])
htmltext+='''<p><b>Recommended structure:</b> one stable section Terrei backbone, plus species-specific extensions or membership views using the same permanent locus identifiers. Keep the section profile for cross-species context and species profiles for fine-scale comparisons. Do not compare raw allele counts across different locus panels without reporting their common callable denominator.</p><p><b>Initial species pilots:</b> citrinoterreus (8 labels), hortae (6) and pseudoterreus (6) have concordant group boundaries in the current core and type-derived three-gene assignments. They are the most defensible first pilots, but their small panels are insufficient for a general species scheme. Terreus/montoensis require the taxonomy/backbone conflict gate first. The one- and two-genome species panels cannot establish population-wide core loci. The 95% and 100% screens coincide for n&lt;20 because a single absence already exceeds 5%.</p><p><b>Fungal allele definition:</b> prefer a pilot based on validated, exon-joined genomic CDS, with one representative isoform per gene. Assess whole genomic loci including introns only as a separately defined allele system. Protein orthogroups identify candidates; they do not define nucleotide alleles. Raw transcript abundance or absence must not be treated as genome-locus absence. chewBBACA documents CDS input for schema creation and allele calling; this offers a testable route after fungal gene-model validation, not a claim that bacterial gene prediction is appropriate for Aspergillus.</p>'''
htmltext+=table(GATES,[('Gate',lambda r:r['gate']),('Status',lambda r:r['status']),('Inputs',lambda r:r['input']),('Work',lambda r:r['action']),('Outputs / acceptance',lambda r:r['output']+' '+r['criterion'])])
htmltext+='''<h2>Next topology analyses, in priority order</h2><ol><li><b>Support-aware core branches:</b> calculate core-locus gene and site concordance factors (gCF/sCF), their discordant alternatives and decisive denominators. First verify locus alignments and orthology. The seven overlapping marker combinations must not be counted as independent gene votes.</li><li><b>Focused hypothesis tests:</b> infer unconstrained and species-monophyly-constrained alternatives for the terreus/montoensis region, then evaluate candidates with an AU test on the same alignment/model. Do this separately for marker and validated core data. Use branch-length optimization and multiple constrained searches; failure to reject is not proof of the constraint.</li><li><b>Low-information and provenance sensitivity:</b> collapse unsupported/near-zero branches, analyze one representative per sequence identity group, vary taxon sampling, and add verified type genomes. Repeat the clade and quartet summaries to see which conflicts persist.</li><li><b>Gene-wise evolutionary conflict:</b> after per-locus QC, compare concatenation with a coalescent species-tree approach and examine localized discordance. Consider recombination/introgression only after excluding paralogy, contamination, sample mix-ups and weak signal; a tanglegram alone cannot diagnose the cause.</li><li><b>Independent validation:</b> retain holdout genomes and evaluate species assignments, cgMLST callability and within-species discrimination without tuning on the test panel.</li></ol><p>These extensions are prepared as scientific gates, not launched as new large jobs. This assessment used one CPU and a 2 GiB hard memory limit for each sequential worker; BUSCO/Panaroo settings were not changed.</p>'''
htmltext+='<h2>Methods and limitations</h2><ul>'+''.join('<li>'+esc(s)+'</li>' for s in MAN['limitations'])+'</ul>'
htmltext+='''<p>Unrooted RF is |S1 symmetric-difference S2| / (|S1|+|S2|), excluding terminal splits. Within-set RF uses induced trees on the same target members. Nearest split maximizes F1 overlap, with precision used as a tie break; extra/missing tips are descriptive candidates, not a uniquely identified set of misclassified genomes. Quartet topologies are determined from four-point sums of unit-edge path distances, with unresolved quartets counted separately. All quartets are resolved in the stored binary point trees; this does not imply strong support. They are non-independent and no significance p-values are assigned. The short-edge RF threshold (1.1e-6) is an explicit numerical sensitivity, not a biological cutoff.</p><p>Trees are display-rooted on LOCAL_SRR10507650. Rooting and rotations preserve all unrooted splits and all 2,145 pairwise patristic distances per tree within 1e-9. The left core order is fixed. Marker child rotations minimize connector inversions conditional on this fixed order. Branch-length axes share a numerical extent within each pair, not across all pages. No connector-crossing statistic is used as biological evidence.</p>'''
htmltext+='<h2>Evidence and primary documentation</h2><ul><li><a href="https://iqtree.github.io/doc/Concordance-Factor">IQ-TREE: gCF and sCF definitions, decisive denominators and computation</a></li><li><a href="https://iqtree.github.io/doc/Advanced-Tutorial#tree-topology-tests">IQ-TREE: topology tests and constrained searches</a></li><li><a href="https://chewbbaca.readthedocs.io/en/latest/user/getting_started/quick_start.html">chewBBACA: CDS input, allele calling and locus-subset selection</a></li><li><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC6857279/">OrthoFinder: phylogenetic orthology and duplication inference</a></li><li><a href="https://pubmlst.org/organisms/aspergillus-fumigatus">PubMLST Aspergillus fumigatus: conventional MLST context and reported limitations</a></li></ul><p>A focused search did not establish an existing curated section Terrei cgMLST scheme. This is not proof that none exists; registry and literature checks remain part of the prerelease nomenclature gate.</p>'
htmltext+='<h2>Download figures</h2>'+gallery
htmltext+='<h2>Machine-readable audits</h2><ul>'+''.join(f'<li><a href="data/{p.name}">{p.name}</a></li>' for p in sorted((OUT/'data').glob('*.json')))+'</ul></body></html>'
(OUT/'review.html').write_text(htmltext,encoding='utf-8')
(OUT/'qa/figure_checks.json').write_text(json.dumps(QA,indent=2));assert all(not q['out_of_page_text'] for q in QA),QA
(OUT/'status.json').write_text(json.dumps({'state':'rendered_pending_visual_review','figures':OUTPUTS,'limits':LIMITS},indent=2))
shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
print('COMPLETE',len(OUTPUTS),flush=True)
