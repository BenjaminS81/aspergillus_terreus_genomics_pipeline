"""Additional plain genome and gene-architecture figures; preserve prior outputs."""
from circos_review_common import *
LIMITS=constrain_process()
import sys, math, csv, hashlib, shutil, gc
os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib_cache')
sys.path.insert(0,str(OUT/'vendor'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch, FancyArrow
from pycirclize import Circos
import pymupdf

P=OUT/'plain_views';P.mkdir(exist_ok=True);(P/'qa').mkdir(exist_ok=True)
D=OUT/'data';order=readjson(D/'figure_order.json');ann=readjson(D/'sample_annotations.json')
metrics=readjson(D/'genome_metrics.json')
by_sample={s:m for m in metrics.values() for s in m['samples']}
taxa=sorted({ann['LOCAL_'+s]['taxon'] for s in order})
taxcolors=dict(zip(taxa,['#8367a0','#b09a32','#477ca6','#8a5747','#c77796','#378e88','#8d969f']))
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'text.color':'#243d50','axes.labelcolor':'#243d50'})
def label(s):
    a=ann['LOCAL_'+s];flag='*' if a['species_ring'] in ('Unresolved','Name conflict') else ''
    return f'{order.index(s)+1:02d}  {a["taxon"]}{flag} | {s}'
def save(fig,stem):
    for ext in ('pdf','svg','png'):fig.savefig(P/(stem+'.'+ext),dpi=110,bbox_inches='tight',pad_inches=.25,facecolor='white')
    plt.close(fig);gc.collect();print('CREATED',stem,flush=True)
def bar(fig,x,cmap,lo,hi,title,ticks,log=False):
    ax=fig.add_axes([x,.115,.14,.009]);norm=LogNorm(lo,hi) if log else Normalize(lo,hi)
    fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=ax,orientation='horizontal',ticks=ticks)
    ax.set_title(title,fontsize=9,pad=5);ax.tick_params(labelsize=8,length=2)

def genome_view():
    c=Circos({s:by_sample[s]['assembly_bp'] for s in order},start=10,end=350,space=.4)
    tracks=[('GC (%)','gc_percent',(87,94),'viridis',40,60),('GC skew','gc_skew',(77,84),'PuOr_r',-.08,.08),('N bases (%)','n_percent',(67,74),'Greys',0,5),('Contig starts','contig_starts',(57,64),'YlOrRd',0,4)]
    for sec in c.sectors:
        s=sec.name;m=by_sample[s]
        sec.text(label(s),r=113,orientation='vertical',size=8)
        sec.add_track((98,100)).axis(fc=taxcolors[ann['LOCAL_'+s]['taxon']],ec='none')
        for name,key,radius,cmap,lo,hi in tracks:
            tr=sec.add_track(radius);cm=plt.get_cmap(cmap);norm=Normalize(lo,hi)
            for b in m['bins']:
                value=b[key]
                if key=='contig_starts':value=math.log10(1+value)
                tr.rect(b['start'],b['end'],fc='#dfe3e6' if value is None else cm(norm(value)),ec='none',lw=0)
        sec.add_track((47,54)).axis(fc=plt.get_cmap('Blues')(LogNorm(1e3,1e7)(m['N50_bp'])),ec='none')
        sec.add_track((45,45.4)).xticks([0,m['assembly_bp']/2,m['assembly_bp']],labels=['','',''],tick_length=.7)
    for name,r in [('Reported taxon',99),('GC content',90.5),('GC skew',80.5),('N bases',70.5),('Contig starts',60.5),('N50',50.5)]:c.text(name,r=r,deg=0,size=8)
    c.text('WHOLE-GENOME OVERVIEW',r=18,deg=0,size=20,weight='bold')
    c.text('77 assembly labels\n72 distinct assemblies\n\nSequence composition\nand assembly continuity',r=3,deg=0,size=15,linespacing=1.7)
    c.text('Sector length = assembly size\n100 positional bins per assembly',r=23,deg=180,size=11,linespacing=1.7)
    fig=plt.figure(figsize=(25,26));ax=fig.add_axes([.075,.22,.85,.66],projection='polar');c.plotfig(ax=ax);ax.set_ylim(0,127)
    fig.text(.055,.968,'Aspergillus | whole-genome Circos',fontsize=27,weight='bold',va='top')
    fig.text(.055,.938,'All 77 canonical assembly labels in one plot; angular length follows total assembly size.\nContigs are ordered by decreasing length within each assembly. This display order does not imply chromosome adjacency or synteny.',fontsize=11,va='top',linespacing=1.6,color='#546b7b')
    fig.text(.075,.157,'OUTSIDE TO INSIDE   Reported taxon  /  GC content  /  GC skew  /  N bases  /  contig-start density  /  N50',fontsize=12,weight='bold')
    for x,args in zip([.075,.255,.435,.615,.795],[('viridis',40,60,'GC (%)',[40,50,60]),('PuOr_r',-.08,.08,'GC skew, clipped',[-.08,0,.08]),('Greys',0,5,'N bases (%)',[0,2.5,5]),('YlOrRd',0,4,'log10(1 + contig starts)',[0,2,4]),('Blues',1e3,1e7,'N50 (bp; log scale)',[1e3,1e5,1e7],True)]):bar(fig,x,*args)
    fig.legend(handles=[Patch(color=taxcolors[t],label=t) for t in taxa],loc='center',bbox_to_anchor=(.5,.068),ncol=7,frameon=False,fontsize=9)
    fig.text(.055,.025,'GC uses A+C+G+T as denominator; N content uses all bases. Skew uses submitted contig orientation and is descriptive. N50 is an assembly-level summary.\nNumerical color limits can saturate; exact metrics and contig offsets are provided. * Reported name remains unresolved or conflicting. Identical-assembly aliases are retained.',fontsize=9,color='#546b7b')
    save(fig,'all_genomes_plain_circos')

def architecture_view():
    data={d['sample']:d for d in readjson(D/'visualization_data.json')}
    panel=readjson(D/'reference_genes.json');genes=[g['gene'] for g in panel if g['role']=='pathway_panel']
    colors=dict(zip(genes,['#5377aa','#157f87','#86b87a','#b5a042','#d48948','#9a70ad','#d36780','#54a5bf','#34675b']))
    rows=[(s,j,p) for s in order for j,p in enumerate(data[s]['segments'])]
    fig,ax=plt.subplots(figsize=(31,51));fig.subplots_adjust(left=.29,right=.98,top=.962,bottom=.047)
    def xy(p,start,end,strand):
        if p['orientation']==1:return start+p['reference_offset'],end+p['reference_offset']+1,strand
        return -end+p['reference_offset'],-start+p['reference_offset']+1,'+' if strand=='-' else '-'
    def arrow(p,g,y,color,height=.17):
        lo=max(p['window_start'],g['start']);hi=min(p['window_end'],g['end'])
        if hi<lo:return None
        a,b,st=xy(p,lo,hi,g['strand']);a/=1000;b/=1000
        clipped=(g['strand']=='+' and g['end']>p['window_end']) or (g['strand']=='-' and g['start']<p['window_start'])
        ax.add_patch(FancyArrow(a if st=='+' else b,y,b-a if st=='+' else a-b,0,width=height,head_width=height*1.25,head_length=0 if clipped else min((b-a)*.4,.4),length_includes_head=True,fc=color,ec='white',lw=.2))
        return (a+b)/2
    for i,(s,j,p) in enumerate(rows):
        y=i*1.5
        if order.index(s)%2==0:ax.axhspan(y-.66,y+.7,color='#f5f7f8',zorder=-5)
        lo,hi,_=xy(p,p['window_start'],p['window_end'],'+');ax.plot([lo/1000,hi/1000],[y,y],color='#a9b9c4',lw=.8)
        for edge in (lo,hi):ax.plot([edge/1000]*2,[y-.12,y+.12],c='#7a8e9c',lw=.6)
        ax.text(-.016,y-.1,label(s) if j==0 else f'    {s} | separate piece {j+1}',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=9,weight='bold' if j==0 else 'normal')
        ax.text(-.016,y+.25,f'{p["contig"]}: {p["window_start"]:,}-{p["window_end"]:,}  ({"+" if p["orientation"]==1 else "-"} display)',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7.2,color='#617786')
        for g in p['annotation_genes']:
            if not any(min(g['end'],h['end'])>=max(g['start'],h['start']) for h in p['genes']+p['flanking_reference_genes']):arrow(p,g,y,'#ccd4da',.12)
        labels=[]
        for kind,key in [('flank','flanking_reference_genes'),('panel','genes')]:
            for g in p[key]:
                col=colors[g['gene']] if kind=='panel' else '#657c8e';x=arrow(p,g,y,col)
                if x is not None:labels.append((x,g['gene'],kind,col))
        lasts=[-1e9]*3
        for x,name,kind,col in sorted(labels):
            spacing=2.0 if kind=='panel' else 4.1;free=[k for k in range(3) if x-lasts[k]>spacing];level=free[0] if free else min(range(3),key=lambda k:lasts[k]);lasts[level]=x
            yy=y-.23-level*.19;ax.plot([x,x],[y-.09,yy+.04],c=col,lw=.35)
            ax.text(x,yy,name if kind=='panel' else name.replace('ATEG_',''),ha='center',va='bottom',fontsize=7.2 if kind=='panel' else 5.9,color=col,weight='bold' if kind=='panel' else 'normal')
    ax.set_xlim(15,96);ax.set_ylim(len(rows)*1.5,-1);ax.set_yticks([]);ax.spines[['left','right','top']].set_visible(False)
    ax.set_xticks(range(20,100,5));ax.grid(axis='x',color='#e8edef',lw=.4,zorder=-10)
    ax.set_xlabel('Genomic spacing in the NIH2624 anchor frame (kb); each contig piece is anchored independently',fontsize=12)
    ax.secondary_xaxis('top')
    fig.text(.035,.991,'Lovastatin | gene architecture and adjacent arrays',fontsize=27,weight='bold',va='top')
    fig.text(.035,.979,'77 assembly labels; 101 independently displayed contig pieces. Gene order, strand and genomic spacing are retained, with up to 20 kb of flanking sequence.\nColoured arrows identify panel homologs; dark grey arrows are reference flank homologs; light grey arrows are other validated annotation models.',fontsize=11,va='top',linespacing=1.7,color='#546b7b')
    fig.legend(handles=[Patch(color=colors[g],label=g) for g in genes]+[Patch(color='#657c8e',label='Reference flank'),Patch(color='#ccd4da',label='Other annotated gene')],loc='lower center',bbox_to_anchor=(.53,.025),ncol=11,frameon=False,fontsize=10)
    fig.text(.035,.009,'Coordinates at left are native 1-based inclusive intervals. Flank labels abbreviate ATEG_09957-09974. A clipped directional end has no arrowhead.\nSeparate rows do not imply adjacency; reference anchoring is not a basewise alignment. Unlinked homologs are included for completeness and do not establish a complete cluster.\nFlanking gene architecture is shown where matching annotation or reference-homology evidence exists. Exact identifiers and provenance are supplied in region_features.tsv. * Unresolved/conflicting reported taxon.',fontsize=10,color='#546b7b')
    save(fig,'lovastatin_plain_linear_architecture')

if __name__=='__main__':
    original={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (OUT/'figures').glob('*') if p.is_file()}
    genome_view();architecture_view()
    checks={'all_77_labels':len(order)==77 and set(order)==set(by_sample),'unique_assemblies_72':len(metrics)==72,'prior_plots_unchanged':all(hashlib.sha256((OUT/p).read_bytes()).hexdigest()==h for p,h in original.items())}
    for stem in ('all_genomes_plain_circos','lovastatin_plain_linear_architecture'):
        doc=pymupdf.open(P/(stem+'.pdf'));text=doc[0].get_text();checks[stem+'_one_page']=len(doc)==1;checks[stem+'_long_identifiers_present']=all(s in text for s in order if len(s)>4)
        if stem=='all_genomes_plain_circos':checks['genome_figure_no_locus_specific_content']=not any(w in text.lower() for w in ('lovastatin','lovb','cluster','pathway'))
        doc[0].get_pixmap(matrix=pymupdf.Matrix(.6 if 'circos' in stem else .4,.6 if 'circos' in stem else .4)).save(str(P/'qa'/(stem+'.png')))
        if 'linear' in stem:
            for tag,y in [('top',0),('middle',1700),('bottom',doc[0].rect.height-400)]:doc[0].get_pixmap(matrix=pymupdf.Matrix(1,1),clip=pymupdf.Rect(0,y,doc[0].rect.width,min(y+400,doc[0].rect.height))).save(str(P/'qa'/('linear_'+tag+'.png')))
        doc.close()
    assert all(checks.values()),checks
    writejson(P/'validation.json',{'checks':checks,'resource_limits':LIMITS,'previous_plot_hashes':original})
    shutil.copy2(D/'genome_metrics.json',P/'genome_metrics.json');shutil.copy2(D/'region_features.tsv',P/'region_features.tsv')
    with (P/'assembly_metrics.tsv').open('w',newline='',encoding='utf-8') as h:
        fields=['sample','taxon','assembly_bp','gc_percent','N50_bp','contigs','n_percent','representative','sequence_sha256'];w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader()
        for s in order:w.writerow({'sample':s,'taxon':ann['LOCAL_'+s]['taxon'],**{k:by_sample[s][k] for k in fields[2:]}})
    shutil.copy2(Path(__file__),P/Path(__file__).name)
    (P/'README.txt').write_text('Additional plain views, 2026-11-13. All previous figures preserved.\n\nall_genomes_plain_circos: all 77 labels (72 distinct assemblies), with reported taxon, GC, skew, Ns, contig-start density and N50. No locus-specific tracks or links.\n\nlovastatin_plain_linear_architecture: all 77 labels, 101 separate pieces, panel genes and adjacent flanking gene arrays with order, orientation and spacing. No QC/status columns.\n\nBoth figures are supplied as vector PDF/SVG and PNG. Source code runs from the original workspace and reads the prior validated data; no source FASTA modification or new alignment is needed. Reproduction: python plot_plain_terreus_views.py\n\nGenomic composition uses 100 bins per assembly. Contigs are sorted by descending length, then ID; this does not infer joins. GC=(G+C)/(A+C+G+T), skew=(G-C)/(G+C), N percentage uses all bases. N50 is computed without length filtering or splitting at Ns. Composition color limits may saturate; JSON preserves exact metrics and display offsets.\n\nThe same 2 GiB worker memory cap and one-CPU affinity are enforced. The prior figure hashes are retained in validation.json.\n',encoding='utf-8')
    print('VALIDATED',sum(checks.values()),'checks; all prior plots preserved',flush=True)
