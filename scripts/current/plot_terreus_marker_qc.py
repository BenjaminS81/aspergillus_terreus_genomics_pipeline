from pathlib import Path
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10';A=P/'analysis';F=P/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
rows=json.loads((A/'orthology_crosscheck/gene_models_vs_genomic_loci.json').read_text());rows=sorted(rows,key=lambda d:d['sample'])
fields=['OG0009417_gene_models','OG0000286_gene_models','OG0002455_gene_models','OG0000542_gene_models','OG0003828_gene_models','BenA_genomic_candidates','CaM_genomic_candidates','RPB2_genomic_candidates','BenA_annotation_models_matching_marker','CaM_annotation_models_matching_marker','RPB2_annotation_models_matching_marker']
labels=['Old tubulin\nOG0009417','BenA-associated\nOG0000286','Fusion-associated\nOG0002455','CaM-associated\nOG0000542','RPB2\nOG0003828','BenA\ngenomic','CaM\ngenomic','RPB2\ngenomic','BenA\nannotation','CaM\nannotation','RPB2\nannotation']
data=np.array([[d[k] for k in fields] for d in rows]);fig,ax=plt.subplots(figsize=(13.3,17.3));cmap=ListedColormap(['#eceff3','#438f8b','#efb35e','#c97749','#963d3b']);norm=BoundaryNorm(np.arange(-.5,5,1),cmap.N)
ax.imshow(data,cmap=cmap,norm=norm,aspect='auto');ax.set_yticks(range(len(rows)),[d['sample'] for d in rows],fontsize=8);ax.set_xticks(range(len(fields)),labels,fontsize=8);ax.xaxis.tick_top();ax.tick_params(length=0,pad=5)
for i in range(len(rows)):
    for j in range(len(fields)):ax.text(j,i,str(data[i,j]),ha='center',va='center',fontsize=7,color='white' if data[i,j] in [1,3,4] else '#283345')
ax.axvline(4.5,color='white',lw=5);ax.axvline(7.5,color='white',lw=5)
fig.suptitle('Orthogroup membership and genomic marker recovery',x=.12,y=.994,ha='left',fontsize=17,fontweight='bold',color='#203550')
fig.text(.12,.970,'66 core-tree assemblies · Counts of gene models, strong genomic candidates and matching annotation models',fontsize=10,color='#526273')
fig.text(.12,.018,'A zero annotation count does not imply genomic absence. Orthogroups can contain fragments, related genes or fused models.\nGenomic candidates use reference-locus BLAST criteria; copy-number claims still require gene-model and read-depth review.',fontsize=9,color='#526273')
fig.subplots_adjust(left=.17,right=.99,top=.91,bottom=.055)
for ext in ['svg','pdf','png']:fig.savefig(F/('orthogroups_vs_genomic_loci.'+ext),dpi=170)
plt.close(fig)
calls=json.loads((A/'genomic_markers/locus_calls.json').read_text());samples=sorted({d['sample'] for d in calls});lookup={(d['sample'],d['locus']):d for d in calls}
fig,axs=plt.subplots(1,3,figsize=(12,17.2),gridspec_kw={'width_ratios':[1,1,1]})
for j,l in enumerate(['BenA','CaM','RPB2']):
    ax=axs[j];cov=[lookup[s,l].get('query_coverage_percent',0) for s in samples];ident=[lookup[s,l].get('weighted_identity_percent',0) for s in samples]
    col=['#d5a240' if lookup[s,l]['status']=='partial_single_candidate' else '#408e89' if lookup[s,l].get('sequence') else '#b7bfc8' for s in samples]
    ax.barh(range(len(samples)),cov,color=col,height=.76);ax.set_xlim(0,111);ax.set_ylim(len(samples)-.4,-.6);ax.set_title(l,loc='left',fontweight='bold',fontsize=13,color='#203550');ax.set_yticks(range(len(samples)),samples if j==0 else ['']*len(samples),fontsize=7);ax.set_xlabel('Reference-query coverage (%)',fontsize=9);ax.tick_params(axis='y',length=0);ax.grid(axis='x',alpha=.16);ax.set_axisbelow(True)
    for i,s in enumerate(samples):ax.text(102,i,f'{ident[i]:.1f}' if ident[i] else 'hold',va='center',fontsize=6.4,color='#526273')
fig.suptitle('Marker recovery across all assembly labels',x=.13,y=.994,ha='left',fontsize=17,fontweight='bold',color='#203550')
fig.text(.13,.970,'77 labels · Teal: near-complete candidate · Amber: partial candidate · Grey: unresolved\nNumbers beside bars show BLAST identity (%). An unresolved screen is not evidence of gene loss.',fontsize=9,color='#526273')
fig.subplots_adjust(left=.18,right=.99,top=.93,bottom=.045,wspace=.18)
for ext in ['svg','pdf','png']:fig.savefig(F/('assembly_marker_recovery.'+ext),dpi=170)
plt.close(fig)
print('Saved marker QC figures',F)
