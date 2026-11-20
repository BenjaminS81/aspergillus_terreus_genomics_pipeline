from pathlib import Path
import json,math,collections,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from terreus_tree_helpers import *
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10';A=P/'analysis';F=P/'figures';F.mkdir(exist_ok=True)
meta={d['reference_id']:d for d in json.loads((A/'tree_reference_metadata.json').read_text())};occ={d['id']:d for d in json.loads((A/'tree_inputs/locus_occupancy.json').read_text())};sources={d['sample']:d for d in json.loads((A/'genome_sources.json').read_text())}
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
def label(ident):
    if ident.startswith('LOCAL_'):
        s=ident[6:];return s+(' [1 locus]' if occ[ident]['available_loci']==1 else '')
    d=meta[ident];s=d['analysis_species'].replace('Aspergillus_','A. ').replace('_',' ')+' '+d['strain']
    if d['type_flag']:s+=' [T]'
    if d.get('historical_type_name'):s+=' [type of '+d['historical_type_name'].replace('Aspergillus_','A. ')+']'
    if ident.startswith('HIST2011'):s+=' (2011)'
    return s
def color(ident):
    if ident.startswith('LOCAL_'):return '#167d78' if sources[ident[6:]]['in_core_tree'] else '#4d99b8'
    d=meta[ident]
    return '#b07b26' if d['role']=='outgroup_reference' else '#8c919b' if d['role']=='historical_reference' else '#223e61' if d['type_flag'] else '#536173'
def positions(root,cladogram=False):
    tips=root.leaves();ys={id(n):i for i,n in enumerate(tips)};xs={id(root):0}
    def visit(n):
        for c in n.children:xs[id(c)]=xs[id(n)]+(1 if cladogram else c.length);visit(c)
        if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
    visit(root);return xs,ys,tips
tree_names=['combined_all_strains','BenA','CaM','RPB2','type_focused_sensitivity']
for name in tree_names:
    if not (A/'trees'/(name+'.iqtree')).exists():continue
    root=parse_newick((A/'trees'/(name+'.treefile')).read_text());xs,ys,tips=positions(root);maxx=max(xs.values());fig,ax=plt.subplots(figsize=(17,max(8,len(tips)*.17)))
    def draw(n):
        if n.children:
            ax.plot([xs[id(n)]]*2,[min(ys[id(c)] for c in n.children),max(ys[id(c)] for c in n.children)],color='#5c6a77',lw=.5)
            for c in n.children:ax.plot([xs[id(n)],xs[id(c)]],[ys[id(c)]]*2,color='#5c6a77',lw=.5);draw(c)
            if strong(n):ax.plot(xs[id(n)],ys[id(n)],'o',color='#253e59',ms=1.6)
        else:ax.text(xs[id(n)]+maxx*.012,ys[id(n)],label(n.name),va='center',fontsize=6.8,color=color(n.name),fontweight='bold' if n.name.startswith('LOCAL_') else 'normal')
    draw(root);ax.set_ylim(len(tips)+1,-2);ax.set_xlim(-maxx*.02,maxx*1.85);ax.set_yticks([]);ax.set_xlabel('Substitutions per site');ax.spines[['left','right','top']].set_visible(False)
    ax.set_title(('Combined BenA + CaM + RPB2' if name=='combined_all_strains' else name)+' reference phylogeny\n'+f'{len(tips)} tips · ML tree · dots: SH-aLRT ≥80 and UFBoot ≥95 · root: A. flavipes NRRL302',loc='left',fontsize=14,pad=20,color='#203550')
    fig.subplots_adjust(left=.03,right=.99,top=.98,bottom=.025)
    for ext in ['svg','pdf','png']:fig.savefig(F/(name+'_reference_phylogram.'+ext),dpi=140,bbox_inches='tight',pad_inches=.2)
    plt.close(fig)
    if name not in ['combined_all_strains','type_focused_sensitivity']:continue
    keep={k for k,d in meta.items() if d['type_flag'] or d.get('accepted_type_name') or d.get('historical_type_name')}|{n.name for n in tips if n.name.startswith('LOCAL_')}
    focused=root if name=='type_focused_sensitivity' else prune(root,keep);xs,ys,tips=positions(focused,True);R=max(xs.values());n=len(tips)
    # Equal-angle cladogram keeps densely sampled strain labels readable; lengths intentionally omitted.
    theta={id(t):2*math.pi*i/n for i,t in enumerate(tips)}
    def angles(node):
        for c in node.children:angles(c)
        if node.children:theta[id(node)]=sum(theta[id(c)] for c in node.children)/len(node.children)
    angles(focused);fig,ax=plt.subplots(figsize=(19,19),subplot_kw={'projection':'polar'})
    def polar(node):
        r=xs[id(node)]/R
        if node.children:
            low=min(theta[id(c)] for c in node.children);high=max(theta[id(c)] for c in node.children)
            ax.plot([low+(high-low)*i/50 for i in range(51)],[r]*51,color='#677887',lw=.7)
            for c in node.children:ax.plot([theta[id(c)]]*2,[r,xs[id(c)]/R],color='#677887',lw=.7);polar(c)
            if strong(node):ax.scatter(theta[id(node)],r,s=6,c='#253e59',zorder=3)
        else:
            a=theta[id(node)];ax.plot([a,a],[r,1.015],color='#c8d1d9',lw=.45,ls=':');degree=math.degrees(a);left=90<degree<270
            text_label=label(node.name)
            if node.name in meta and meta[node.name].get('historical_type_name'):
                d=meta[node.name];text_label=d['strain']+' [type: '+d['historical_type_name'].replace('Aspergillus_','A. ')+']'
            ax.text(a,1.07,text_label,rotation=degree+(180 if left else 0),rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=7.7,color=color(node.name),fontweight='bold' if node.name.startswith('LOCAL_') else 'normal')
            loci=occ[node.name]['available_loci'];ax.scatter(a,1.035,s=13,c='#438f8b' if loci==3 else '#e7ac50' if loci==2 else '#bc6550',marker='s',zorder=4)
    polar(focused);ax.set_ylim(0,1.72);ax.set_axis_off();fig.suptitle('Terrei reference placement',fontsize=24,x=.04,y=.98,ha='left',color='#203550',fontweight='bold')
    subtitle=f'{len(tips)} tips · Independent type-focused sensitivity analysis' if name=='type_focused_sensitivity' else f'{len(tips)} type-focused tips, pruned from the full {len(root.leaves())}-tip tree'
    fig.text(.04,.949,subtitle+'\nCircular cladogram: branch lengths are not represented. Project labels are retained.',fontsize=11,color='#536173')
    handles=[Line2D([0],[0],color=c,lw=3,label=l) for l,c in [('Core-tree assembly','#167d78'),('Additional assembly label','#4d99b8'),('Type reference [T]','#223e61'),('Outgroup','#b07b26')]]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.035,.025),frameon=False,ncol=2,fontsize=10)
    fig.text(.60,.04,'Outer squares: 3 loci (teal), 2 loci (amber), 1 locus (red)\nDots: SH-aLRT ≥80 and UFBoot ≥95\nTaxonomic placement remains subject to single-locus concordance.',fontsize=9,color='#536173')
    fig.subplots_adjust(top=.92,bottom=.07,left=.02,right=.98)
    basename='combined_type_focused_circular' if name=='type_focused_sensitivity' else 'combined_type_focused_pruned'
    for ext in ['svg','pdf','png']:fig.savefig(F/(basename+'.'+ext),dpi=180,bbox_inches='tight',pad_inches=.25)
    plt.close(fig)
print('Rendered all completed trees')
