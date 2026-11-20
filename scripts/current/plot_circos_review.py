"""Three distinct data-backed Circos proposals and a fully stacked locus comparison."""
from circos_review_common import *
LIMITS=constrain_process()
import sys, math, gc, csv, hashlib, shutil
os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib_cache')
sys.path.insert(0,str(OUT/'vendor'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrow
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from pycirclize import Circos
import pycirclize
D=OUT/'data';F=OUT/'figures'
data={x['sample']:x for x in readjson(D/'visualization_data.json')}
order=readjson(D/'figure_order.json');ann=readjson(D/'sample_annotations.json')
metrics=readjson(D/'genome_metrics.json');regions=readjson(D/'region_metrics.json')
panel=readjson(D/'reference_genes.json');genes=[g['gene'] for g in panel if g['role']=='pathway_panel']
palette=dict(zip(genes,['#5377aa','#157f87','#86b87a','#b5a042','#d48948','#9a70ad','#d36780','#54a5bf','#34675b']))
taxa=sorted({x['taxon'] for x in data.values()});taxcolors=dict(zip(taxa,['#8367a0','#b09a32','#477ca6','#8a5747','#c77796','#378e88','#8d969f']))
statuscolors={'nine_gene_block':'#22877d','partial_or_divergent_block':'#d89a39','no_convincing_colocalized_block':'#929ca5'}
statusnames={'nine_gene_block':'Nine-gene block','partial_or_divergent_block':'Partial / split','no_convincing_colocalized_block':'No convincing block'}
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'font.size':10,'axes.labelcolor':'#243d50','text.color':'#243d50'})
inventory={};coordmaps=[]
def label(s):
    a=ann['LOCAL_'+s];flag='*' if a['species_ring'] in ('Unresolved','Name conflict') else ''
    typ=' [T]' if a.get('genome_type_species') else ''
    return f'{order.index(s)+1:02d}  {data[s]["taxon"]}{flag} | {s}{typ}'
def save(fig,stem,samples,segments=None):
    for ext in ('pdf','svg','png'):
        fig.savefig(F/(stem+'.'+ext),dpi=110,bbox_inches='tight',pad_inches=.25,facecolor='white')
    inventory[stem]={'samples':list(samples),'sample_count':len(samples),'segment_count':segments,'formats':['pdf','svg','png']}
    plt.close(fig);gc.collect();print('FIGURE',stem,flush=True)
def header(fig,title,subtitle):
    fig.text(.055,.968,title,fontsize=26,weight='bold',va='top')
    fig.text(.055,.938,subtitle,fontsize=11,va='top',linespacing=1.6,color='#546b7b')
def basecircos(sectors):
    c=Circos(sectors,start=10,end=350,space=.4)
    return c
def render(c):
    fig=plt.figure(figsize=(25,26));ax=fig.add_axes([.075,.22,.85,.66],projection='polar')
    c.plotfig(ax=ax);ax.set_ylim(0,127)
    return fig,ax
def labels_and_rings(c):
    for sector in c.sectors:
        s=sector.name;d=data[s]
        sector.text(label(s),r=113,orientation='vertical',size=8)
        sector.add_track((99,100)).axis(fc=taxcolors[d['taxon']],ec='none')
        sector.add_track((96,98)).axis(fc=statuscolors[d['block_status']],ec='none')
def categorical_legend(fig,y=.075):
    handles=[Patch(color=taxcolors[t],label=t) for t in taxa]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,y),ncol=7,frameon=False,fontsize=9)
    fig.legend(handles=[Patch(color=statuscolors[k],label=statusnames[k]) for k in statuscolors],loc='center',bbox_to_anchor=(.5,y-.021),ncol=3,frameon=False,fontsize=10)
def cb(fig,rect,cmap,vmin,vmax,title,ticks=None,log=False):
    ax=fig.add_axes(rect);norm=LogNorm(vmin,vmax) if log else Normalize(vmin,vmax)
    bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=ax,orientation='horizontal',ticks=ticks)
    ax.set_title(title,fontsize=9,pad=5);ax.tick_params(labelsize=8,length=2)
def heat(track,rows,key,cmap,vmin,vmax,offset=0,reverse=False,length=None):
    cm=plt.get_cmap(cmap);norm=Normalize(vmin,vmax)
    for b in rows:
        val=b[key]
        if val is None:color='#dfe3e6'
        else:
            if reverse and key=='gc_skew':val=-val
            color=cm(norm(val))
        x,y=b['start'],b['end']
        if reverse:x,y=length-y,length-x
        track.rect(offset+x,offset+y,fc=color,ec='none',lw=0)
def windowxy(p,start,end):
    if p['orientation']==1:return start-p['window_start'],end-p['window_start']+1
    return p['window_end']-end,p['window_end']-start+1

def cluster_circos():
    sizes={s:sum(p['window_end']-p['window_start']+1 for p in data[s]['segments'])+2000*(len(data[s]['segments'])-1) for s in order}
    c=basecircos(sizes);labels_and_rings(c)
    for sector in c.sectors:
        s=sector.name;d=data[s];tracks={key:sector.add_track(rad) for key,rad in [('plus',(90,94)),('minus',(85,89)),('cov',(80,83)),('dna',(76,79)),('aa',(72,75)),('gc',(62,69)),('skew',(53,60))]}
        offset=0
        for j,p in enumerate(d['segments']):
            r=regions[d['representative']][j];length=r['bp']
            for tr in tracks.values():tr.rect(offset,offset+length,fc='#f0f2f3',ec='none')
            for kind,key in [('annotation','annotation_genes'),('flank','flanking_reference_genes'),('panel','genes')]:
                for g in p[key]:
                    lo=max(p['window_start'],g['start']);hi=min(p['window_end'],g['end'])
                    if lo>hi:continue
                    a,b=windowxy(p,lo,hi);strand=g['strand'] if p['orientation']==1 else ('-' if g['strand']=='+' else '+')
                    tr=tracks['plus' if strand=='+' else 'minus'];color=palette[g['gene']] if kind=='panel' else '#617486' if kind=='flank' else '#cbd2d8'
                    x,y=(offset+a,offset+b) if strand=='+' else (offset+b,offset+a)
                    cut_head=(g['strand']=='+' and g['end']>p['window_end']) or (g['strand']=='-' and g['start']<p['window_start'])
                    tr.arrow(x,y,head_length=0 if cut_head else .16,shaft_ratio=.7,fc=color,ec='white',lw=.08)
                    if kind=='panel':
                        tracks['cov'].rect(offset+a,offset+b,fc=plt.get_cmap('Blues')(Normalize(.7,1)(g['query_coverage'])),ec='none')
                        k='dna' if g['method']=='blastn' else 'aa';vmin=75 if k=='dna' else 45
                        tracks[k].rect(offset+a,offset+b,fc=plt.get_cmap('magma')(Normalize(vmin,100)(g['identity_percent'])),ec='none')
            heat(tracks['gc'],r['bins'],'gc_percent','viridis',30,70,offset,p['orientation']==-1,length)
            heat(tracks['skew'],r['bins'],'gc_skew','PuOr_r',-.5,.5,offset,p['orientation']==-1,length)
            coordmaps.append({'sample':s,'segment':j+1,'contig':p['contig'],'genomic_start_1based':p['window_start'],'genomic_end_1based_inclusive':p['window_end'],'display_start_0based':offset,'display_end_0based_exclusive':offset+length,'orientation':p['orientation']})
            offset+=length+2000
    c.text('A | LOCUS CONTEXT',r=27,deg=0,size=23,weight='bold')
    c.text('77 assembly labels\n101 separate contig segments\nUp to 20 kb of flanking sequence',r=16,deg=0,size=14,linespacing=1.6)
    c.text('Outside to inside\nReported taxon / block status\nForward / reverse gene models\nQuery coverage / DNA identity / AA identity\nGC content / oriented GC skew\n\nWhite radial gaps separate contig pieces.\nGrey feature tracks = annotated / flank homologs.\nBackground alone does not establish absence.',r=12,deg=180,size=11,linespacing=1.7)
    fig,ax=render(c);header(fig,'A  |  Lovastatin regions and flanking context','Every assembly in one coordinate-based plot. Angular length follows extracted sequence length; separate pieces have a 2 kb display gap.\nDirection is oriented to the NIH2624 anchor per fragment. The detected-gene span is not a validated biosynthetic-cluster boundary.')
    fig.legend(handles=[Patch(color=palette[g],label=g) for g in genes],loc='center',bbox_to_anchor=(.5,.145),ncol=9,frameon=False,fontsize=11)
    cb(fig,[.075,.107,.14,.009],'Blues',70,100,'Query coverage (%)',[70,85,100])
    cb(fig,[.255,.107,.14,.009],'magma',75,100,'DNA identity (%)',[75,90,100])
    cb(fig,[.435,.107,.14,.009],'magma',45,100,'Protein identity (%)',[45,70,100])
    cb(fig,[.615,.107,.14,.009],'viridis',30,70,'GC, 1 kb (%)',[30,50,70])
    cb(fig,[.795,.107,.13,.009],'PuOr_r',-.5,.5,'GC skew, clipped',[-.5,0,.5])
    categorical_legend(fig,.068)
    fig.text(.055,.017,'* Reported name remains unresolved or conflicting in the marker review. [T] Documented genome type anchor. Homology does not demonstrate production.\nGC = (G+C)/(A+C+G+T). Skew = (G-C)/(G+C); sign changes with displayed orientation. Values outside color limits saturate; exact values accompany the plots.',fontsize=9,color='#546b7b')
    save(fig,'A_lovastatin_regions_circos',order,101)

def genomepos(s,g):
    m=metrics[data[s]['representative']];cc=next(x for x in m['contig_coordinates'] if x['contig']==g['contig'])
    return cc['start_0based']+g['start']-1,cc['start_0based']+g['end']
def whole_genome_circos():
    c=basecircos({s:metrics[data[s]['representative']]['assembly_bp'] for s in order});labels_and_rings(c)
    for sector in c.sectors:
        s=sector.name;d=data[s];m=metrics[d['representative']]
        tr=sector.add_track((89,94));heat(tr,m['bins'],'gc_percent','viridis',40,60)
        tr=sector.add_track((81,87));heat(tr,m['bins'],'gc_skew','PuOr_r',-.08,.08)
        tr=sector.add_track((73,79));heat(tr,m['bins'],'n_percent','Greys',0,5)
        tr=sector.add_track((65,71));tr.heatmap([np.log10(1+b['contig_starts']) for b in m['bins']],vmin=0,vmax=4,cmap='YlOrRd',rect_kws={'ec':'none','lw':0})
        tr=sector.add_track((58,63));tr.axis(fc='#eef1f3',ec='none')
        for p in d['segments']:
            for g in p['genes']:
                x,y=genomepos(s,g);tr.rect(x,y,fc=palette[g['gene']],ec='none')
                # A visible tick shows the midpoint; it does not enlarge the interval.
                tr.scatter([(x+y)/2],[.5],vmin=0,vmax=1,s=3,c=palette[g['gene']],marker='|')
        tr=sector.add_track((49,55));tr.axis(fc=plt.get_cmap('Blues')(LogNorm(1000,10000000)(m['N50_bp'])),ec='none')
        # Relative sector positions are assembly display coordinates, never chromosome positions.
        sector.add_track((47.5,48)).xticks([0,m['assembly_bp']/2,m['assembly_bp']],labels=['','',''],tick_length=.8)
    ref='NIH2624';rhit=next(g for p in data[ref]['segments'] for g in p['genes'] if g['gene']=='lovB');ra,rb=genomepos(ref,rhit)
    links=0
    for s in order:
        if s==ref or data[s]['block_status']!='nine_gene_block':continue
        hit=next(g for p in data[s]['segments'] for g in p['genes'] if g['gene']=='lovB');a,b=genomepos(s,hit)
        c.link_line((ref,(ra+rb)/2),(s,(a+b)/2),r1=46,r2=46,color=palette['lovB'],alpha=.45,lw=.65);links+=1
    fig,ax=render(c);header(fig,'B  |  All assemblies in genomic coordinates','77 sectors; angular length is proportional to assembly size. Contigs are ordered by decreasing length within each assembly.\nThis is an assembly display coordinate system: adjacent contigs are not inferred to be joined or syntenic. Each composition track has 100 bins per assembly.')
    ax.text(0,0,'77 assemblies\n72 distinct sequence sets\n\nLinks: lovB homologs\nin nine-gene blocks',ha='center',va='center',fontsize=14,bbox={'facecolor':'white','alpha':.93,'edgecolor':'none','pad':10})
    fig.text(.08,.156,'OUTSIDE → INSIDE   Taxon · block status · GC · GC skew · N bases · contig-start density · pathway loci · N50',fontsize=12,weight='bold')
    cb(fig,[.075,.112,.14,.009],'viridis',40,60,'GC (%)',[40,50,60])
    cb(fig,[.255,.112,.14,.009],'PuOr_r',-.08,.08,'GC skew, clipped',[-.08,0,.08])
    cb(fig,[.435,.112,.14,.009],'Greys',0,5,'N bases (%)',[0,2.5,5])
    cb(fig,[.615,.112,.14,.009],'YlOrRd',0,4,'log10(1 + contig starts)',[0,2,4])
    cb(fig,[.795,.112,.13,.009],'Blues',1000,10000000,'N50 (bp; log scale)',[1e3,1e5,1e7],True)
    categorical_legend(fig,.07)
    fig.text(.055,.018,f'{links} links connect NIH2624 lovB to detected homologs in other nine-gene-block labels; these are locus-homology links, not genome-wide alignment ribbons.\nSkew uses submitted contig orientation and is descriptive only. Color limits may saturate; tables retain exact values. * Unresolved/conflicting reported name. [T] Documented type anchor.',fontsize=9,color='#546b7b')
    save(fig,'B_all_genomes_coordinate_circos',order)
    inventory['B_all_genomes_coordinate_circos']['lovB_links']=links

def dashboard_circos():
    c=basecircos({s:1 for s in order});labels_and_rings(c)
    specs=[('Size (Mb)','assembly_bp',1e6,24,36,'viridis'),('GC (%)','gc_percent',1,49,54,'viridis'),('N50 (log bp)','N50_bp',1,1e3,1e7,'Blues'),('Contigs (log)','contigs',1,1,15000,'YlOrRd'),('N bases (%)','n_percent',1,0,5,'Greys')]
    for sector in c.sectors:
        s=sector.name;d=data[s];m=metrics[d['representative']]
        for j,(name,key,scale,lo,hi,cm) in enumerate(specs):
            norm=LogNorm(lo,hi) if key in ('N50_bp','contigs') else Normalize(lo,hi)
            sector.add_track((89-j*6,94-j*6)).axis(fc=plt.get_cmap(cm)(norm(m[key]/scale)),ec='none')
        found={g['gene']:g for p in d['segments'] for g in p['genes']}
        for j,gene in enumerate(genes):
            hit=found.get(gene);color='#e4e8eb'
            if hit:
                rgb=np.array(matplotlib.colors.to_rgb(palette[gene]));cov=min(1,hit['query_coverage']);color=tuple(1-(1-rgb)*cov)
            tr=sector.add_track((59-j*4,62-j*4));tr.axis(fc=color,ec='none')
            if hit and hit['method']=='tblastn':tr.scatter([.5],[.5],vmin=0,vmax=1,marker='.',s=7,c='#172e3e')
        typ=bool(ann['LOCAL_'+s].get('genome_type_species'));sector.add_track((20,22)).axis(fc='#222e46' if typ else '#eef1f3',ec='none')
    for j,(name,*_) in enumerate(specs):c.text(name,r=91.5-j*6,deg=0,size=8)
    for j,g in enumerate(genes):c.text(g,r=60.5-j*4,deg=0,size=9,weight='bold',color=palette[g])
    c.text('Type anchor',r=21,deg=0,size=7)
    c.text('C | COMPARISON\n77 assembly labels',r=5,deg=180,size=16,weight='bold',linespacing=1.6)
    fig,ax=render(c);header(fig,'C  |  Assembly quality and lovastatin evidence','Equal sectors make every genome equally visible. Read clockwise: the existing 66-genome three-marker order, followed by 11 additional labels.\nFive numerical QC tracks, nine locus-evidence tracks, block status, reported taxon and documented type anchors; this plot contains no inferred tree.')
    fig.text(.08,.153,'Gene colour = strong homolog; shade encodes query coverage (70–100%). Grey = below threshold / not detected. Dot = protein-rescue hit.',fontsize=12)
    for j,(name,key,scale,lo,hi,cm) in enumerate(specs):
        ticks=[lo,hi] if key not in ('N50_bp','contigs') else ([1e3,1e5,1e7] if key=='N50_bp' else [10,100,1000,10000])
        cb(fig,[.075+j*.18,.11,.14,.009],cm,lo,hi,name,ticks,key in ('N50_bp','contigs'))
    categorical_legend(fig,.069)
    fig.text(.055,.017,'Recommended overview for comparing samples. The linear companion resolves individual gene order and coordinates. No convincing block does not prove absence.\nAliases of identical assemblies remain visible and are identified in the accompanying table. * Unresolved/conflicting reported name. [T] Documented genome type anchor.',fontsize=9,color='#546b7b')
    save(fig,'C_all_genomes_qc_evidence_circos',order)

def displayxy(p,start,end,strand):
    if p['orientation']==1:return start+p['reference_offset'],end+p['reference_offset']+1,strand
    return -end+p['reference_offset'],-start+p['reference_offset']+1,'+' if strand=='-' else '-'
def linear_arrow(ax,a,b,strand,y,color,height=.16,complete_head=True):
    a/=1000;b/=1000
    ax.add_patch(FancyArrow(a if strand=='+' else b,y,(b-a) if strand=='+' else (a-b),0,width=height,head_width=height*1.35,head_length=min((b-a)*.4,.35) if complete_head else 0,length_includes_head=True,fc=color,ec='white',lw=.25))
def stacked_linear():
    rows=[(s,j,p) for s in order for j,p in enumerate(data[s]['segments'])]
    fig,ax=plt.subplots(figsize=(34,55));fig.subplots_adjust(left=.27,right=.84,top=.965,bottom=.045)
    for i,(s,j,p) in enumerate(rows):
        d=data[s];y=i*1.55
        if order.index(s)%2==0:ax.axhspan(y-.65,y+.87,color='#f4f7f8',zorder=-5)
        lo,hi,_=displayxy(p,p['window_start'],p['window_end'],'+');corelo,corehi,_=displayxy(p,p['start'],p['end'],'+')
        ax.plot([lo/1000,hi/1000],[y,y],color='#bac6ce',lw=1)
        ax.plot([corelo/1000,corehi/1000],[y,y],color='#dce5e9',lw=8,zorder=-1,solid_capstyle='butt')
        for edge in (lo,hi):ax.plot([edge/1000]*2,[y-.15,y+.15],color='#627786',lw=.8)
        left=label(s) if j==0 else f'    {s} | piece {j+1}'
        ax.text(-.016,y-.12,left,transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=9,weight='bold' if j==0 else 'normal')
        coordinate=f'{p["contig"]}: {p["window_start"]:,}-{p["window_end"]:,}  ({"+" if p["orientation"]==1 else "-"} display)'
        ax.text(-.016,y+.26,coordinate,transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7.4,color='#617786')
        for g in p['annotation_genes']:
            if any(min(g['end'],h['end'])>=max(g['start'],h['start']) for h in p['genes']+p['flanking_reference_genes']):continue
            a,b,st=displayxy(p,max(g['start'],p['window_start']),min(g['end'],p['window_end']),g['strand']);cut=(g['strand']=='+' and g['end']>p['window_end']) or (g['strand']=='-' and g['start']<p['window_start']);linear_arrow(ax,a,b,st,y,'#cbd3d9',.1,not cut)
        labelled=[]
        for kind,key in [('flank','flanking_reference_genes'),('panel','genes')]:
            for g in p[key]:
                a,b,st=displayxy(p,max(g['start'],p['window_start']),min(g['end'],p['window_end']),g['strand'])
                if b<=a:continue
                color=palette[g['gene']] if kind=='panel' else '#617486';cut=(g['strand']=='+' and g['end']>p['window_end']) or (g['strand']=='-' and g['start']<p['window_start']);linear_arrow(ax,a,b,st,y,color,complete_head=not cut)
                labelled.append(((a+b)/2000,g['gene'],kind,color))
        # Place names on alternating levels, with explicit leader lines to their features.
        lasts=[-1e9,-1e9,-1e9]
        for x,name,kind,color in sorted(labelled):
            width=2.1 if kind=='panel' else 4.3
            free=[k for k in range(3) if x-lasts[k]>width];level=free[0] if free else min(range(3),key=lambda k:lasts[k]);lasts[level]=x
            yy=y-.24-level*.19;ax.plot([x,x],[y-.1,yy+.05],color=color,lw=.4)
            ax.text(x,yy,name if kind=='panel' else name.replace('ATEG_',''),ha='center',va='bottom',fontsize=7.1 if kind=='panel' else 5.8,color=color,weight='bold' if kind=='panel' else 'normal')
        if j==0:
            right=f'{statusnames[d["block_status"]]}\n{d["genome_wide_strong_genes"]}/9 strong genome-wide; {d["segment_count"]} piece(s)'
            ax.text(1.014,y,right,transform=ax.get_yaxis_transform(),va='center',fontsize=8,color=statuscolors[d['block_status']],linespacing=1.6)
        else:ax.text(1.014,y,'Independent contig segment',transform=ax.get_yaxis_transform(),va='center',fontsize=7.5,color='#617786')
        if p['fragment_edge']:ax.text(hi/1000+.25,y,'!',va='center',fontsize=11,color='#bf6546',weight='bold')
    ax.set_ylim(len(rows)*1.55,-1);ax.set_xlim(15,96);ax.set_yticks([]);ax.spines[['left','right','top']].set_visible(False)
    ax.set_xticks(range(20,100,5));ax.grid(axis='x',color='#e5eaed',lw=.5,zorder=-10)
    ax.set_xlabel('NIH2624 anchor frame (kb); every contig piece is positioned independently',fontsize=13)
    ax.secondary_xaxis('top').set_xlabel('Gene arrows preserve genomic span and direction; grey background marks the span of detected panel genes',fontsize=11)
    fig.text(.035,.991,'Lovastatin loci | all 77 assemblies with flanking regions',fontsize=28,weight='bold',va='top')
    fig.text(.035,.981,'101 contig segments; up to 20 kb flanks on each side, clipped at contig ends. Each panel homolog and detected reference flank is labelled.\nCoordinates at left are native 1-based inclusive assembly coordinates. Flank labels 09957–09974 abbreviate ATEG_ identifiers; grey arrows are other annotation models.',fontsize=12,va='top',color='#546b7b')
    fig.legend(handles=[Patch(color=palette[g],label=g) for g in genes]+[Patch(color='#617486',label='Reference flank homolog'),Patch(color='#cbd3d9',label='Other annotation model')],loc='lower center',bbox_to_anchor=(.53,.023),ncol=11,frameon=False,fontsize=11)
    fig.text(.035,.009,'! A detected panel gene lies within 2 kb of a contig edge. White space between rows never implies genomic adjacency. Reference anchoring is not a basewise alignment.\nThe span of detected genes does not define a complete BGC. Scattered homologs and negative screens do not establish pathway function or biological absence.\nAll model identifiers, exact coordinates, identity, query coverage and provenance are in data/region_features.tsv. * Unresolved/conflicting reported taxon. [T] Documented genome type anchor.',fontsize=11,color='#546b7b')
    save(fig,'lovastatin_all_genomes_stacked_flanks',order,len(rows))

def key():
    fig,axes=plt.subplots(1,2,figsize=(16,17))
    for ax in axes:ax.axis('off')
    for i,s in enumerate(order):
        col=i//39;row=i%39;d=data[s];m=metrics[d['representative']]
        text=label(s)+('\n       Identical sequence aliases: '+', '.join(m['samples']) if len(m['samples'])>1 else '')
        axes[col].text(0,1-row*.0246,text,fontsize=7.6,va='top',linespacing=1.3)
    fig.suptitle('Genome key | shared order in all four comparisons',fontsize=21,x=.055,ha='left')
    fig.subplots_adjust(left=.055,right=.99,top=.94,bottom=.025)
    save(fig,'genome_key',order)

if __name__=='__main__':
    actions={'A':cluster_circos,'B':whole_genome_circos,'C':dashboard_circos,'linear':stacked_linear,'key':key}
    selected=sys.argv[1:] or list(actions)
    old=readjson(D/'figure_inventory.json') if (D/'figure_inventory.json').exists() else {}
    for name in selected:actions[name]()
    old.update(inventory);writejson(D/'figure_inventory.json',old)
    if coordmaps:writejson(D/'region_display_coordinates.json',coordmaps)
    writejson(D/'plotting_environment.json',{'python':sys.version,'pycirclize':pycirclize.__version__,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'resource_limits':LIMITS})
