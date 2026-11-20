from pathlib import Path
import json,math,collections,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow,Patch,Polygon
from matplotlib.backends.backend_pdf import PdfPages
from terreus_tree_helpers import parse_newick
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';L=P/'lovastatin';F=P/'figures/lovastatin';F.mkdir(parents=True,exist_ok=True);IND=F/'individual';IND.mkdir(exist_ok=True)
data={x['sample']:x for x in json.loads((L/'visualization_data.json').read_text())};panel=json.loads((L/'reference_genes.json').read_text());ref={x['gene']:x for x in panel};genes=[x['gene'] for x in panel if x['role']=='pathway_panel'];palette=dict(zip(genes,['#5377aa','#157f87','#86b87a','#b5a042','#d48948','#9a70ad','#d36780','#54a5bf','#34675b']))
root=parse_newick((P/'core66_combinations/BenA_CaM_RPB2.treefile').read_text());order=[t.name[6:] for t in root.leaves()];order+=sorted(set(data)-set(order));assert len(order)==77
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
def export(fig,stem,pdf=True):
    for ext in (['svg','pdf','png'] if pdf else ['svg','png']):fig.savefig(F/(stem+'.'+ext),dpi=155,bbox_inches='tight',pad_inches=.2)
    plt.close(fig)
def displayxy(p,c):
    if p['orientation']==1:return c['start']+p['reference_offset'],c['end']+p['reference_offset'],c['strand']
    return -c['end']+p['reference_offset'],-c['start']+p['reference_offset'],'+' if c['strand']=='-' else '-'
def arrow(ax,start,end,strand,y,color,height=.44):
    start/=1000;end/=1000;w=end-start
    ax.add_patch(FancyArrow(start if strand=='+' else end,y,w if strand=='+' else -w,0,width=height,head_width=height,head_length=min(abs(w)*.4,.65),length_includes_head=True,facecolor=color,edgecolor='white',lw=.25,zorder=3))
def overview(flanks):
    tracks=[]
    for s in order:
        for j,p in enumerate(data[s]['segments'] or [None]):tracks.append((s,j,p))
    fig,ax=plt.subplots(figsize=(25,max(16,len(tracks)*.23+4)));allx=[]
    for i,(s,j,p) in enumerate(tracks):
        label=(data[s]['taxon']+' | '+s) if j==0 else '  contig fragment '+str(j+1)
        ax.text(-.018,i,label,transform=ax.get_yaxis_transform(),va='center',ha='right',fontsize=7.8,color='#213b54',weight='bold' if j==0 else 'normal')
        if p is None:ax.text(36,i,'No strong candidate under the stated screen',va='center',fontsize=7,color='#8c9299');continue
        lo,hi=(p['window_start'],p['window_end']) if flanks else (p['start'],p['end']);a,b,_=displayxy(p,{'start':lo,'end':hi,'strand':'+'});allx.extend([a/1000,b/1000]);ax.plot([a/1000,b/1000],[i,i],c='#b8c4cd',lw=.8)
        for edge in [a/1000,b/1000]:ax.plot([edge]*2,[i-.19,i+.19],c='#687a88',lw=.6)
        if flanks:
            for c in p['annotation_genes']:
                if any(min(c['end'],g['end'])>=max(c['start'],g['start']) for g in p['genes']):continue
                x1,x2,st=displayxy(p,c);arrow(ax,x1,x2,st,i,'#d5dbe0',.3)
            for c in p['flanking_reference_genes']:
                x1,x2,st=displayxy(p,c);arrow(ax,x1,x2,st,i,'#8595a3',.34)
        for c in p['genes']:
            x1,x2,st=displayxy(p,c);arrow(ax,x1,x2,st,i,palette[c['gene']])
        if p['fragment_edge']:ax.text(b/1000+.6,i,'edge',va='center',fontsize=6,color='#b86a43')
    ax.set_ylim(len(tracks),-2);ax.set_xlim(min(allx)-2,max(allx)+7);ax.set_yticks([]);ax.spines[['left','right','top']].set_visible(False);ax.set_xlabel('kb in the NIH2624 reference frame; each contig fragment anchored independently')
    ax.set_title('Lovastatin-region homologs | '+('with up to 20 kb flanks' if flanks else 'nine-gene pathway panel')+'\n77 assembly labels; separate rows preserve contig breaks and unlinked candidates',loc='left',fontsize=19,weight='bold',color='#213b54',pad=25)
    fig.legend(handles=[Patch(color=palette[g],label=g) for g in genes]+([Patch(color='#8595a3',label='Reference flank homolog'),Patch(color='#d5dbe0',label='Other annotated gene')] if flanks else []),loc='lower center',ncol=6,frameon=False,fontsize=10)
    fig.subplots_adjust(left=.28,right=.98,top=.95,bottom=.065);fig.text(.28,.032,'Fragments are not joined. Scattered homologs do not establish a biosynthetic cluster. Arrow lengths use genomic spans; function and expression are untested.',fontsize=9,color='#536778');export(fig,'all_genomes_linear'+('_20kb_flanks' if flanks else '_core'))
overview(False);overview(True)
# Circular overview is an explicit evidence matrix, not a circular chromosome claim.
fig,ax=plt.subplots(figsize=(22,22),subplot_kw={'projection':'polar'});start=.16;span=2*math.pi-.65;width=span/9;inner=18
for i,s in enumerate(order):
    d=data[s];found={g['gene']:g for p in d['segments'] for g in p['genes']}
    for j,g in enumerate(genes):
        c=found.get(g);col=palette[g] if c and d['block_status']=='nine_gene_block' else '#d49a62' if c else '#e4e7ea'
        ax.bar(start+(j+.5)*width,.87,width=width*.97,bottom=inner+i,color=col,edgecolor='white',lw=.15)
    ax.text(start-.02,inner+i+.4,f'{i+1:02d}',ha='right',va='center',fontsize=5.5,rotation=math.degrees(start)-90,rotation_mode='anchor',color='#30465a')
for j,g in enumerate(genes):a=start+(j+.5)*width;ax.text(a,inner+80,g,ha='center',va='center',fontsize=15,weight='bold',color=palette[g])
ax.set_ylim(0,inner+87);ax.axis('off');ax.text(0,0,'77 assembly labels\n9 pathway genes',ha='center',va='center',fontsize=16,color='#213b54');fig.suptitle('Lovastatin pathway evidence across all genomes',x=.045,y=.975,ha='left',fontsize=25,weight='bold',color='#213b54')
fig.text(.045,.94,'Circular evidence matrix: numbered tracks match the accompanying 77-genome key; sectors are genes.\nGene colour = homolog in a contiguous nine-gene block; amber = homolog in a partial/dispersed context; grey = below threshold/not detected.',fontsize=11,color='#536778');fig.text(.045,.025,'A missing hit is not proof of biological absence. Matching sequence does not demonstrate lovastatin production. Repeated assembly aliases are retained and documented.',fontsize=10,color='#536778');export(fig,'all_genomes_circular_evidence')
fig,axes=plt.subplots(1,2,figsize=(14,14))
for ax in axes:ax.axis('off')
for i,s in enumerate(order):
    col=i//39;row=i%39;d=data[s];axes[col].text(0,1-row*.024,f'{i+1:02d}  {d["taxon"]} | {s}',fontsize=9,color='#213b54',va='top')
fig.suptitle('Circular lovastatin comparison | genome key',x=.055,y=.98,ha='left',fontsize=20,weight='bold',color='#213b54');fig.subplots_adjust(top=.94,bottom=.02,left=.055,right=.98);export(fig,'circular_genome_key')
if '--overview-only' in sys.argv:sys.exit(0)
def linear_detail(ax,d,flanks):
    pieces=d['segments']
    if not pieces:ax.text(.5,.5,'No strong candidate',transform=ax.transAxes,ha='center');ax.axis('off');return
    for j,p in enumerate(pieces):
        lo,hi=(p['window_start'],p['window_end']) if flanks else (p['start'],p['end']);base=lo
        ax.plot([0,(hi-lo)/1000],[j,j],c='#899ba8',lw=.9);ax.text(-.025,j,p['contig'],transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7)
        if flanks:
            for c in p['annotation_genes']:arrow(ax,c['start']-base,c['end']-base,c['strand'],j,'#d9dfe3',.22)
            for c in p['flanking_reference_genes']:arrow(ax,c['start']-base,c['end']-base,c['strand'],j,'#8f9da9',.28)
        previous=-1e9;level=0
        for c in sorted(p['genes'],key=lambda x:x['start']):
            mid=(c['start']+c['end']-2*base)/2000;level=(1-level) if mid-previous<3 else 0;previous=mid
            arrow(ax,c['start']-base,c['end']-base,c['strand'],j,palette[c['gene']],.35);ax.text(mid,j-.28-.13*level,c['gene'],fontsize=7,ha='center')
        ax.text((hi-lo)/1000+.5,j,f'{lo:,}-{hi:,}',fontsize=6.5,va='center',color='#647788')
    mx=max(((p['window_end']-p['window_start']) if flanks else p['end']-p['start'])/1000 for p in pieces);ax.set_xlim(-.5,mx+max(8,mx*.22));ax.set_ylim(len(pieces)-.4,-.65);ax.set_yticks([]);ax.set_xlabel('kb from the displayed segment start; native genomic orientation',fontsize=8);ax.spines[['left','right','top']].set_visible(False);ax.tick_params(labelsize=8)
def circ_detail(ax,d,flanks):
    pieces=d['segments'];radius=1
    if not pieces:ax.text(0,0,'No strong candidate',ha='center');ax.axis('off');return
    lengths=[(p['window_end']-p['window_start']+1) if flanks else (p['end']-p['start']+1) for p in pieces];total=sum(lengths);gap=.15;usable=2*math.pi-gap*len(pieces);offset=.1
    for i,(p,length) in enumerate(zip(pieces,lengths)):
        base=p['window_start'] if flanks else p['start'];arc=usable*length/total;factor=arc/length
        ax.bar(offset+arc/2,.02,width=arc,bottom=radius,color='#8b9da9',lw=0)
        others=([({**c,'gene':'annotation'},'#d5dce2') for c in p['annotation_genes']]+[(c,'#8d9ba8') for c in p['flanking_reference_genes']]) if flanks else []
        lastlabel=-1e9;labellevel=0
        for c,col in others+[(c,palette[c['gene']]) for c in sorted(p['genes'],key=lambda x:x['start'])]:
            lo=max(0,c['start']-base);hi=min(length,c['end']-base+1)
            if hi<=lo:continue
            a=offset+lo*factor;b=offset+hi*factor;ax.bar((a+b)/2,.11,width=b-a,bottom=.95,color=col,edgecolor='white',lw=.3)
            end=b if c['strand']=='+' else a;back=end-min(.035,(b-a)*.35) if c['strand']=='+' else end+min(.035,(b-a)*.35)
            ax.add_patch(Polygon([[back,.95],[end,1.005],[back,1.06]],closed=True,color=col,zorder=5))
            if c.get('gene') in palette:
                mid=(a+b)/2;labellevel=(1-labellevel) if mid-lastlabel<.18 else 0;lastlabel=mid
                ax.text(mid,1.13+.09*labellevel,c['gene'],fontsize=7,ha='center',va='center',color=col)
        ax.text(offset+arc/2,.8,'segment '+str(i+1),ha='center',fontsize=7,color='#50687a');offset+=arc+gap
    ax.set_ylim(0,1.35);ax.axis('off');ax.text(0,0,'Linear contig segments\nwrapped for display\nGaps = separate segments',ha='center',va='center',fontsize=8,color='#536778')
atlas=F/'all_genomes_cluster_and_flanks_atlas.pdf'
with PdfPages(atlas) as pdf:
    for i,s in enumerate(order):
        d=data[s];fig=plt.figure(figsize=(17,13));grid=fig.add_gridspec(2,2,height_ratios=[.8,1.2]);a=fig.add_subplot(grid[0,0]);b=fig.add_subplot(grid[0,1]);c=fig.add_subplot(grid[1,0],projection='polar');e=fig.add_subplot(grid[1,1],projection='polar')
        linear_detail(a,d,False);linear_detail(b,d,True);circ_detail(c,d,False);circ_detail(e,d,True)
        for ax,title in [(a,'Core homologs | linear'),(b,'Core + up to 20 kb surroundings | linear'),(c,'Core homologs | circular display'),(e,'Core + surroundings | circular display')]:ax.set_title(title,loc='left',fontsize=12,color='#213b54',pad=12)
        fig.suptitle(d['taxon']+' | '+s,fontsize=22,x=.045,y=.98,ha='left',weight='bold',color='#213b54');fig.text(.045,.945,f"{d['block_status'].replace('_',' ')} | {d['genome_wide_strong_genes']}/9 strong genome-wide homologs | {d['segment_count']} displayed segments",fontsize=11,color='#536778')
        fig.legend(handles=[Patch(color=palette[g],label=g) for g in genes]+[Patch(color='#8f9da9',label='Reference flank'),Patch(color='#d9dfe3',label='Other annotated gene')],loc='lower center',bbox_to_anchor=(.5,.045),ncol=6,frameon=False,fontsize=9)
        fig.text(.045,.018,f'Page {i+1}/77. Homology evidence only; scattered homologs are not a cluster. Circular maps do not imply a circular chromosome. Coordinates and evidence: gene_evidence.tsv.',fontsize=8,color='#536778')
        fig.subplots_adjust(top=.89,bottom=.13,left=.12,right=.96,wspace=.3,hspace=.25);pdf.savefig(fig,bbox_inches='tight')
        for ext in ['svg','png']:fig.savefig(IND/(s+'.'+ext),dpi=110,bbox_inches='tight')
        plt.close(fig)
        if i%10==0:print('LOV FIGURE',i+1,'/77',s,flush=True)
(L/'figure_order.json').write_text(json.dumps(order,indent=2));print('LOVASTATIN ATLAS COMPLETE',flush=True)
