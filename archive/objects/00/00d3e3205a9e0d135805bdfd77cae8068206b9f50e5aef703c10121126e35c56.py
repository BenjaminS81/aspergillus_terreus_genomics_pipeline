"""Focused, version-labelled core versus combined-marker comparisons."""
from pathlib import Path
import ast,os,json,hashlib,shutil,gc,html,itertools,sys
ROOT=Path(__file__).resolve().parent
from circos_review_common import constrain_process
LIMITS=constrain_process()
from terreus_tree_helpers import Node,parse_newick,prune
import numpy as np
OUT=ROOT/'stakeholder_release/2026-11-20_v5_core_vs_combined_versions'
OLD=ROOT/'stakeholder_release/2026-11-20_v4_clade_congruence_cgmlst'
for d in ('figures','data','qa','scripts'):(OUT/d).mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(OUT/'qa/mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import fitz
ANN={a['tip']:a for a in json.loads((OLD/'data/core_tip_assignments.json').read_text())};NAMES=sorted(ANN);ALL=set(NAMES)
PAL=json.loads((ROOT/'stakeholder_release/2026-11-17_v3_taxonomy_review/data/palette.json').read_text())['species']
INK='#203849';GRAY='#70818f';GOLD='#d9a72f';OUTPUTS=[];QA=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,'pdf.fonttype':42,'svg.fonttype':'none'})
# Reuse only audited function definitions, never the prior scripts' top-level jobs.
functions={
 'extend_clade_congruence.py':{'nodes','canon','splitmap','rf','distances','reroot','layout_order','newick'},
 'render_clade_congruence.py':{'positions','branchplot','output'}}
for fn,names in functions.items():
 tree=ast.parse((ROOT/fn).read_text(encoding='utf-8'));selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
 assert {n.name for n in selected}==names
 exec(compile(ast.Module(body=selected,type_ignores=[]),fn,'exec'),globals())
CORE=OLD/'data/Core.matched66.treefile';core=parse_newick(CORE.read_text());ORDER=[t.name for t in core.leaves()];rank={t:i for i,t in enumerate(ORDER)}
SPECS=[
 ('01_core_vs_combined_66_de_novo',ROOT/'comparative_visualization/2026-11-11/core66_combinations/BenA_CaM_RPB2.treefile','66-genome de novo analysis','All three loci; inferred directly on the 66-genome cohort.','07_tanglegram_core_vs_BenA_CaM_RPB2'),
 ('02_core_vs_combined_230_complete',ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene/runs/complete230.treefile','Current complete 230-tip analysis','75 genomes + 155 references with all three loci; current taxonomic reference context.','08_tanglegram_core_vs_Three_gene_reference230'),
 ('03_core_vs_combined_283_historical',ROOT/'taxonomy_reference_build/2026-11-10/analysis/trees/combined_all_strains.treefile','Historical 283-tip analysis','76 genomes + 207 references; includes sparse-locus tips. Historical comparison, not the current complete analysis.',None),
 ('04_core_vs_combined_108_type_focused',ROOT/'taxonomy_reference_build/2026-11-10/analysis/trees/type_focused_sensitivity.treefile','108-tip type-focused sensitivity analysis','75 genomes + 33 references; altered reference/outgroup sampling, with fixed partition models from the full historical analysis. Sensitivity comparison.',None)]
audit=[]
for stem,path,label,scope,existing in SPECS:
 raw=parse_newick(path.read_text());source_names=[t.name for t in raw.leaves()];assert len(source_names)==len(set(source_names)) and ALL<=set(source_names)
 projected=prune(raw,ALL);projected.length=0;before=set(splitmap(projected));d=distances(projected)
 marker=reroot(projected);layout_order(marker,rank)
 assert set(splitmap(marker))==before and set(t.name for t in marker.leaves())==ALL
 delta=float(np.max(np.abs(d-distances(marker))));assert delta<1e-9
 (OUT/'data'/(stem+'.matched66.treefile')).write_text(newick(marker)+';')
 groups={s:{t for t in ALL if ANN[t]['assigned_species']==s} for s in sorted({a['assigned_species'] for a in ANN.values()})}
 rec={'figure':stem,'version':label,'scope':scope,'source_path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_tips':len(source_names),'source_genome_tips':sum(t.startswith('LOCAL_') for t in source_names),'matched_tips':66,'excluded_core_tips':[],'normalized_RF_to_core':rf(core,marker),'unrooted_topology_preserved':True,'max_patristic_change':delta,'group_recovery':{s:None if len(m)<2 else canon(m,ALL) in splitmap(marker) for s,m in groups.items()}}
 if existing:
  for ext in ('pdf','svg','png'):shutil.copy2(OLD/'figures'/(existing+'.'+ext),OUT/'figures'/(stem+'.'+ext))
  shutil.copy2(OLD/'qa'/(existing+'.png'),OUT/'qa'/(stem+'.png'));OUTPUTS.append(stem)
  prior=next(q for q in json.loads((OLD/'qa/figure_checks.json').read_text()) if q['figure']==existing);prior['figure']=stem;QA.append(prior)
  source_matched=parse_newick((OLD/'data'/('BenA_CaM_RPB2.matched66.treefile' if source_names.__len__()==66 else 'Three_gene_reference230.matched66.treefile')).read_text())
  assert set(splitmap(marker))==set(splitmap(source_matched))
 else:
  lt,lx,ly=positions(core);rt,rx,ry=positions(marker);extent=max(max(lx.values()),max(rx.values()))*1.04
  fig=plt.figure(figsize=(30,19));left=fig.add_axes([.025,.065,.155,.86]);right=fig.add_axes([.82,.065,.155,.86]);labL=fig.add_axes([.187,.065,.246,.86]);links=fig.add_axes([.435,.065,.13,.86]);labR=fig.add_axes([.57,.065,.245,.86])
  branchplot(left,core,lx,ly,extent);branchplot(right,marker,rx,ry,extent);right.invert_xaxis()
  for ax in (labL,labR,links):ax.set_ylim(65.65,-.65);ax.set_xlim(0,1);ax.axis('off')
  yL={t.name:ly[id(t)] for t in lt};yR={t.name:ry[id(t)] for t in rt}
  for t in ORDER:
   a=ANN[t];co=PAL[a['assigned_species']];links.plot([0,1],[yL[t],yR[t]],color=co,alpha=.65,lw=.85)
   for ax,y in ((labL,yL[t]),(labR,yR[t])):
    ax.add_patch(Rectangle((.001,y-.34),.009,.68,color=co,lw=0));ax.text(.032,y,a['display_label'],fontsize=8.4,va='center',weight='bold' if a['is_type'] else 'normal')
    if a['is_type']:ax.scatter(.02,y,marker='*',s=42,c=GOLD,edgecolors='#71551b',lw=.4)
  fig.text(.025,.974,'Core genome vs BenA + CaM + RPB2 | '+label,fontsize=23,weight='bold')
  fig.text(.025,.947,f'{len(source_names)}-tip source reduced to 66 identical genome labels | fixed core order | shared branch-length scale',fontsize=12,color=GRAY)
  left.set_title('CORE GENOME',fontsize=13,loc='left',pad=17);right.set_title('COMBINED THREE-GENE TREE',fontsize=13,loc='right',pad=17)
  output(fig,stem)
 audit.append(rec)
for ext in ('pdf','svg','png'):shutil.copy2(OLD/'figures'/('12_shared_legend.'+ext),OUT/'figures'/('00_shared_legend.'+ext))
shutil.copy2(CORE,OUT/'data/Core.matched66.treefile');shutil.copy2(OLD/'data/core_tip_assignments.json',OUT/'data/core_tip_assignments.json')
checks=[]
for rec in audit:
 with fitz.open(OUT/'figures'/(rec['figure']+'.pdf')) as pdf:
  assert len(pdf)==1;lines=pdf[0].get_text().splitlines();assert all(sum(line.strip()==a['display_label'] for line in lines)==2 for a in ANN.values())
  assert sum('TYPE:' in line for line in lines)==6
  checks.append({'figure':rec['figure'],'full_labels_verified':132,'type_labels_verified':6,'pages':1})
 assert hashlib.sha256(Path(rec['source_path']).read_bytes()).hexdigest()==rec['sha256']
assert all(not q['out_of_page_text'] for q in QA)
(OUT/'data/version_comparison.json').write_text(json.dumps(audit,indent=2));(OUT/'qa/checks.json').write_text(json.dumps({'figures':QA,'labels':checks,'limits':LIMITS},indent=2))
body='''<!doctype html><html><head><meta charset="utf-8"><title>Core vs combined genes: all four versions</title><style>body{font:16px/1.5 system-ui;color:#203849;max-width:1350px;margin:40px auto;padding:20px}table{border-collapse:collapse;width:100%}td,th{border:1px solid #cbd5dd;padding:12px;text-align:left}img{max-width:100%;display:block;margin:20px auto}article{padding:20px;border:1px solid #cbd5dd;margin:25px 0}a{color:#126a8b}</style></head><body><h1>Core vs combined BenA + CaM + RPB2: four versions</h1><p>Every comparison uses exactly the same 66 genome labels, including all three documented core type genomes. Source tree tips are removed only to obtain this identical cohort; no trees are re-inferred. The core topology, original taxon labels and fixed current species-color mapping are unchanged.</p><p>The prior package contained versions 1 and 2. This focused collection adds historical version 3 and sensitivity version 4. These are sampling/inference alternatives, not four independent biological replicates. Colors carry the current three-gene assignments and do not claim that every historical tree independently supports them.</p><table><tr><th>Version</th><th>Source tips</th><th>Compared tips</th><th>Normalized RF</th><th>Scope</th></tr>'''
for r in audit:body+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['version'],r['source_tips'],66,f"{r['normalized_RF_to_core']:.3f}",r['scope']])+'</tr>'
body+='</table><p>RF measures unrooted split disagreement; the current core has no comparable bootstrap evidence. Branch lengths are substitutions per site, not time; the floccosus tip is a display root only. Each pair shares a numerical branch-length scale.</p><p><a href="figures/00_shared_legend.pdf">Shared legend (one PDF)</a> | <a href="data/version_comparison.json">Source hashes, scope and clade recovery audit</a></p>'
for r in audit:
 s=r['figure'];body+=f'<article><h2>{html.escape(r["version"])}</h2><p>{html.escape(r["scope"])}</p><a href="figures/{s}.pdf">Single-page vector PDF</a> | <a href="figures/{s}.png">300-dpi PNG</a> | <a href="figures/{s}.svg">SVG</a><img src="qa/{s}.png"></article>'
body+='</body></html>';(OUT/'review.html').write_text(body,encoding='utf-8')
shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
for fn in functions:shutil.copy2(ROOT/fn,OUT/'scripts'/fn)
(OUT/'status.json').write_text(json.dumps({'state':'pending_visual_review','comparisons':4,'new_comparisons':2,'same_strains':66,'limits':LIMITS},indent=2))
print(json.dumps([{'version':r['version'],'n':r['matched_tips'],'RF':r['normalized_RF_to_core']} for r in audit],indent=2))
