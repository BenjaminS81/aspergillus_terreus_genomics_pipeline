"""Provenance-labelled views of the seven existing reference-inclusive trees."""
from circos_review_common import *
LIMITS=constrain_process()
import sys, math, csv, hashlib, shutil, collections, gc, html, zipfile
os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib_cache')
sys.path.insert(0,str(OUT/'vendor'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from terreus_tree_helpers import parse_newick,strong
import pymupdf

P=OUT/'marker_provenance';F=P/'figures';D=P/'data';Q=P/'qa'
for path in (F,D,Q):path.mkdir(parents=True,exist_ok=True)
A=OLD/'analysis';T=BASE/'taxonomy_reference_build/2026-11-10/analysis'
allann=readjson(A/'tree_annotations_all_combinations.json');allclades=readjson(A/'type_clades_all_combinations.json')
refs={x['reference_id']:x for x in readjson(T/'tree_reference_metadata.json')}
genomes={'LOCAL_'+x['sample']:x for x in readjson(T/'genome_sources.json')}
colors=readjson(A/'taxon_palette.json');sourcecolors={'genome':'#087f8c','reference':'#bd642d'}
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'text.color':'#213c50','font.size':9})
names=['BenA','CaM','RPB2','BenA_CaM','BenA_RPB2','CaM_RPB2','combined_all_strains']
titles={x:x.replace('_',' + ') for x in names};titles['combined_all_strains']='BenA + CaM + RPB2'
annotations={name:{x['tip']:x for x in allann if x['tree']==name} for name in names}
rows=[];summaries=[];filehashes={}

def provenance(k):
    if k in genomes:return 'genome'
    if k in refs:return 'reference'
    raise ValueError('Unclassified sequence '+k)
def label(k,a):
    origin=provenance(k);prefix='[G]' if origin=='genome' else '[R11]' if k.startswith('HIST2011_') else '[R]'
    # Preserve the reviewed taxon, strain, unresolved-name and type flags verbatim.
    old=a['display_label'].replace(' [2011]','')
    return prefix+' '+old
def ringcols(a):
    return [sourcecolors[provenance(a['tip'])],{1:'#d7794c',2:'#e6ad51',3:'#408f87'}.get(a['available_loci'],'#eeeeee'),colors[a['species_ring']],colors[a['type_clade_category']],colors.get(a.get('genome_type_clade_category','No supported type clade'),'#f3f4f6')]
def positions(root,clado=False):
    tips=root.leaves();ys={id(t):i for i,t in enumerate(tips)};xs={id(root):0}
    def walk(n):
        for c in n.children:xs[id(c)]=xs[id(n)]+(1 if clado else c.length);walk(c)
        if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
    walk(root);return xs,ys,tips
def legend(fig,y=.025):
    handles=[Line2D([],[],marker='s',color=sourcecolors['genome'],lw=0,markersize=8,label='[G] Extracted from a genome'),Line2D([],[],marker='o',color=sourcecolors['reference'],lw=0,markersize=8,label='[R] Curated taxonomic reference'),Line2D([],[],marker='o',color=sourcecolors['reference'],lw=0,markersize=8,label='[R11] Historical 2011 study')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.03,y),ncol=3,frameon=False,fontsize=12)
def save(fig,stem,pdf):
    pdf.savefig(fig,bbox_inches='tight',pad_inches=.3)
    for ext in ('svg','png'):fig.savefig(F/(stem+'.'+ext),dpi=100,bbox_inches='tight',pad_inches=.3)
    plt.close(fig);gc.collect()
def draw_linear(root,name,ann,pdf):
    xs,ys,tips=positions(root);mx=max(xs.values());n=len(tips)
    fig=plt.figure(figsize=(27,max(12,n*.18+4)))
    ax=fig.add_axes([.035,.065,.35,.86]);lab=fig.add_axes([.405,.065,.405,.86],sharey=ax);tracks=fig.add_axes([.825,.065,.15,.86],sharey=ax)
    def draw(node):
        if node.children:
            ax.plot([xs[id(node)]]*2,[min(ys[id(c)] for c in node.children),max(ys[id(c)] for c in node.children)],c='#788795',lw=.6)
            for c in node.children:ax.plot([xs[id(node)],xs[id(c)]],[ys[id(c)]]*2,c='#788795',lw=.6);draw(c)
            if strong(node):ax.scatter(xs[id(node)],ys[id(node)],s=4,c='#213c50')
        else:
            origin=provenance(node.name);y=ys[id(node)];ax.plot([xs[id(node)],mx*1.012],[y,y],c='#cbd4dc',lw=.4,ls=':')
            ax.scatter(mx*1.022,y,marker='s' if origin=='genome' else 'o',s=13,c=sourcecolors[origin],zorder=4)
    draw(root)
    for i,t in enumerate(tips):
        a=ann[t.name];origin=provenance(t.name)
        lab.text(0,i,label(t.name,a),va='center',fontsize=8,color=sourcecolors[origin],weight='bold' if origin=='genome' else 'normal')
        for j,col in enumerate(ringcols(a)):tracks.add_patch(plt.Rectangle((j-.38,i-.47),.76,.94,color=col,lw=0))
        tracks.text(4.65,i,a['type_clade'],va='center',fontsize=7,color='#30455a')
    for axis in (ax,lab,tracks):axis.set_ylim(n+.5,-2)
    ax.set_xlim(0,mx*1.05);ax.set_yticks([]);ax.spines[['top','right','left']].set_visible(False);ax.set_xlabel('Substitutions per site (uncompressed)')
    lab.set_xlim(0,1);lab.axis('off');tracks.set_xlim(-.7,6.8);tracks.set_xticks(range(5),['Source','Loci','Species','Ref. type clade','Genome type'],rotation=45,ha='right');tracks.xaxis.tick_top();tracks.set_yticks([]);tracks.spines[:].set_visible(False)
    counts=collections.Counter(provenance(t.name) for t in tips)
    fig.text(.035,.982,titles[name]+' | sequence provenance',fontsize=25,weight='bold',va='top')
    fig.text(.035,.964,f'{n} tips: {counts["genome"]} genome-derived + {counts["reference"]} taxonomic-reference tips. Original topology and branch lengths retained.\nDots on internal nodes: SH-aLRT >=80 and UFBoot >=95. Leaf symbols and label prefixes identify sequence origin.',fontsize=11,va='top',linespacing=1.6,color='#526778')
    legend(fig,.023)
    fig.text(.035,.012,'[T] Reference ex-type; [TG] documented genome ex-type. Type status is independent of sequence origin. * Unresolved reported taxon; [name flag] source-name conflict.',fontsize=10,color='#526778')
    save(fig,name+'_provenance_linear',pdf)
def draw_circular(root,name,ann,clades,pdf):
    xs,ys,tips=positions(root,True);n=len(tips);mx=max(xs.values());step=2*math.pi/n;theta={id(t):i*step for i,t in enumerate(tips)}
    def calc(node):
        for c in node.children:calc(c)
        if node.children:theta[id(node)]=sum(theta[id(c)] for c in node.children)/len(node.children)
    calc(root)
    fig,ax=plt.subplots(figsize=(35,35),subplot_kw={'projection':'polar'});idx={t.name:i for i,t in enumerate(tips)};painted=set()
    for clade in sorted([x for x in clades if x['status']=='single_type_species'],key=lambda x:-len(x['descendant_tips'])):
        inds=[idx[k] for k in clade['descendant_tips'] if k in idx]
        if len(inds)<2 or set(inds)&painted:continue
        lo,hi=min(inds),max(inds)
        if hi-lo+1!=len(inds):continue
        ax.bar((lo+hi)*step/2,1.02,width=len(inds)*step,bottom=0,color=colors[clade['type_species']],alpha=.055,ec='none');painted.update(inds)
    def draw(node):
        r=xs[id(node)]/mx
        if node.children:
            lo=min(theta[id(c)] for c in node.children);hi=max(theta[id(c)] for c in node.children)
            ax.plot([lo+(hi-lo)*j/50 for j in range(51)],[r]*51,c='#788795',lw=.6)
            for c in node.children:ax.plot([theta[id(c)]]*2,[r,xs[id(c)]/mx],c='#788795',lw=.6);draw(c)
            if strong(node):ax.scatter(theta[id(node)],r,s=4,c='#213c50')
        else:
            a=theta[id(node)];d=ann[node.name];origin=provenance(node.name)
            ax.plot([a,a],[r,1.025],c='#cbd4dc',lw=.45,ls=':');ax.scatter(a,1.025,marker='s' if origin=='genome' else 'o',s=11,c=sourcecolors[origin],zorder=5)
            for j,col in enumerate(ringcols(d)):ax.bar(a,.034,width=step*.93,bottom=1.055+j*.046,color=col,lw=0,zorder=3)
            deg=math.degrees(a);left=90<deg<270
            ax.text(a,1.32,label(node.name,d),rotation=deg+(180 if left else 0),rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=8,color=sourcecolors[origin],weight='bold' if origin=='genome' else 'normal')
    draw(root)
    runs=[]
    for i,t in enumerate(tips):
        code=ann[t.name]['type_clade']
        if runs and runs[-1][2]==code:runs[-1][1]=i
        else:runs.append([i,i,code])
    for lo,hi,code in runs:
        if code=='unresolved' or hi-lo<2:continue
        a=(lo+hi)*step/2;deg=math.degrees(a);ax.text(a,.997,code,fontsize=8,ha='center',va='center',rotation=deg-90+(180 if 90<deg<270 else 0),rotation_mode='anchor')
    ax.set_ylim(0,2.03);ax.axis('off');counts=collections.Counter(provenance(t.name) for t in tips)
    fig.text(.035,.982,titles[name]+' | sequence provenance',fontsize=29,weight='bold',va='top')
    fig.text(.035,.959,f'{n} tips: {counts["genome"]} genome-derived + {counts["reference"]} taxonomic-reference tips.\nCircular cladogram: branch lengths are not represented; the linear companion preserves them.',fontsize=13,va='top',linespacing=1.6,color='#526778')
    fig.text(.035,.071,'RINGS: sequence origin | available loci | reviewed species | supported reference type clade | genome-type evidence\n[T] Reference ex-type; [TG] documented genome ex-type. Type status is separate from origin.\n* Unresolved reported species; [name flag] reference-name conflict. Type-clade codes refer to the tree-specific key.',fontsize=12,color='#526778',linespacing=1.6)
    legend(fig,.025);fig.subplots_adjust(left=.02,right=.98,top=.92,bottom=.11)
    save(fig,name+'_provenance_circular',pdf)

def main():
    previous={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (OLD/'figures').rglob('*') if p.is_file() and p.suffix in ('.pdf','.svg','.png')}
    with PdfPages(P/'marker_trees_sequence_provenance.pdf') as pdf:
        for name in names:
            source=(OLD/'reference_pair_combinations'/(name+'.treefile')) if name.count('_')==1 and name!='combined_all_strains' else OLD/'sources/input_trees'/(name+'.treefile')
            raw=source.read_bytes();filehashes[str(source)]=hashlib.sha256(raw).hexdigest();root=parse_newick(raw.decode());ann=annotations[name];tips=[t.name for t in root.leaves()]
            assert set(tips)==set(ann) and len(tips)==len(set(tips))
            clades=[x for x in allclades if x['tree']==name]
            counts=collections.Counter(provenance(k) for k in tips);counts['historical_2011']=sum(k.startswith('HIST2011_') for k in tips)
            for k in tips:
                a=ann[k];origin=provenance(k);m=refs.get(k,{})
                row={'tree':name,'tip':k,'display_label':label(k,a),'sequence_origin':origin,'source_collection':'canonical genome extraction' if origin=='genome' else '2011 taxonomic study' if k.startswith('HIST2011_') else 'ICPA curated reference-sequence table 2026-10-01','source_paper_doi':m.get('source_paper',''),'reference_table_row':m.get('source_row',''),'assembly_path':genomes.get(k,{}).get('windows_path',''),'reference_ex_type':a['type'],'genome_ex_type':bool(a.get('genome_type_species')),'taxon':a['taxon'],'strain':a['strain'],'reviewed_species':a['species_ring'],'type_clade':a['type_clade'],'type_clade_species':a['type_clade_species']}
                for gene in ('BenA','CaM','RPB2'):row[gene+'_accession']=m.get(gene+'_original_accession') or m.get(gene+'GB','')
                rows.append(row)
            draw_circular(root,name,ann,clades,pdf);draw_linear(root,name,ann,pdf)
            summaries.append({'tree':name,'title':titles[name],'tips':len(tips),**counts,'input_sha256':filehashes[str(source)],'pages':[2*len(summaries)+1,2*len(summaries)+2]})
            shutil.copy2(source,D/(name+'.treefile'));print('RENDERED',name,dict(counts),flush=True)
    writejson(D/'provenance_annotations.json',rows);writejson(D/'tree_summary.json',summaries)
    writejson(D/'type_clade_key.json',allclades);shutil.copy2(A/'taxon_palette.json',D/'taxon_palette.json')
    with (D/'provenance_annotations.tsv').open('w',newline='',encoding='utf-8') as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    checks={'seven_reference_inclusive_trees':len(summaries)==7,'all_tips_classified':all(r['sequence_origin'] in sourcecolors for r in rows),'source_trees_unchanged':all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in filehashes.items()),'prior_figures_unchanged':all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in previous.items()),'genome_ex_types_remain_genome_derived':all(r['sequence_origin']=='genome' for r in rows if r['genome_ex_type']),'reference_ex_types_remain_reference_derived':all(r['sequence_origin']=='reference' for r in rows if r['reference_ex_type'])}
    book=pymupdf.open(P/'marker_trees_sequence_provenance.pdf');checks['fourteen_pdf_pages']=len(book)==14
    for i,page in enumerate(book):
        text=page.get_text();checks[f'page_{i+1}_both_sources_labelled']='[G]' in text and '[R]' in text
        page.get_pixmap(matrix=pymupdf.Matrix(.35,.35)).save(str(Q/f'page_{i+1:02d}.png'))
        if i in (0,1,12,13):
            rect=page.rect;clip=pymupdf.Rect(0,0,rect.width,min(rect.height,650)) if i%2 else pymupdf.Rect(rect.width*.6,rect.height*.28,rect.width*.97,rect.height*.70)
            page.get_pixmap(matrix=pymupdf.Matrix(.85,.85),clip=clip).save(str(Q/f'page_{i+1:02d}_detail.png'))
    book.close();assert all(checks.values()),checks
    writejson(D/'validation.json',{'checks':checks,'resource_limits':LIMITS,'input_tree_hashes':filehashes})
    shutil.copy2(Path(__file__),P/'plot_marker_provenance.py')
    # A compact offline gallery makes every combination accessible without hiding provenance.
    options=''.join(f'<option value="{n}">{html.escape(titles[n])}</option>' for n in names)
    table=''.join('<tr>'+''.join(f'<td>{html.escape(str(s[k]))}</td>' for k in ('title','tips','genome','reference','historical_2011'))+'</tr>' for s in summaries)
    doc='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Marker trees | sequence provenance</title><style>body{font:16px system-ui;color:#213c50;background:#f2f6f8;margin:35px auto;max-width:1250px;padding:0 24px}p{line-height:1.65}a{color:#087f8c}select{padding:10px;font-size:17px;margin:12px}img{display:block;max-width:100%;background:white}table{border-collapse:collapse;background:white}td,th{padding:10px 20px;border-bottom:1px solid #d8e2e8}code{font-weight:bold}.box{background:white;padding:20px;border-radius:10px}.g{color:#087f8c}.r{color:#bd642d}</style><h1>Single genes and combinations | sequence provenance</h1><p>New views of the seven existing reference-inclusive phylogenies. Tree topology, branch lengths, reviewed labels and type-clade calls are retained; previous plots are preserved.</p><div class="box"><b class="g">Square / [G] / teal: genome-extracted sequence.</b><br><b class="r">Circle / [R] / orange: curated taxonomic reference. [R11]: historical 2011 study.</b><p>Provenance is distinct from type status: [T] identifies reference ex-types; [TG] identifies documented ex-type genome cultures. The taxonomic references provide anchors, with unresolved labels and source-name conflicts retained rather than treated as infallible ground truth.</p></div><p><a href="marker_trees_sequence_provenance.pdf">Download the 14-page vector PDF atlas</a> | <a href="data/provenance_annotations.tsv">Tip provenance, accessions and source assemblies</a> | <a href="data/type_clade_key.json">Tree-specific type-clade key</a></p><select id="tree" onchange="showPlot()">'''+options+'''</select><select id="layout" onchange="showPlot()"><option value="circular">Circular</option><option value="linear">Linear</option></select><p><a id="svg" href="figures/BenA_provenance_circular.svg">Open scalable SVG</a> | <a id="png" href="figures/BenA_provenance_circular.png">Open PNG</a></p><img id="figure" src="figures/BenA_provenance_circular.png" alt="Tree with genome and reference sequence provenance"><h2>Coverage</h2><table><tr><th>Tree</th><th>Total tips</th><th>Genome-derived</th><th>Taxonomic references</th><th>2011 subset</th></tr>'''+table+'''</table><p>The curated reference collection is the archived ICPA reference-sequence table dated 2026-10-01, accession-checked in the earlier analysis. The historical collection is explicitly linked to DOI <a href="https://doi.org/10.3114/sim.2011.69.04">10.3114/sim.2011.69.04</a>. The curated table is a collection source, not an assertion that every sequence was newly generated by the 2024 paper.</p><p>Rings/columns: sequence origin; number of available loci; reviewed species; supported reference type clade; genome-type evidence. Type-clade support requires SH-aLRT >=80 and UFBoot >=95. Grey species indicates unresolved assignments; name conflicts remain flagged. Circular trees are cladograms; linear trees retain uncompressed branch lengths.</p><p>Reference-inclusive trees are shown because both origins are present. The prior 66-genome-only comparisons remain available unchanged; all tips in those reduced trees are genome-derived. No sequence search, alignment or tree inference was rerun.</p><script>function showPlot(){const stem='figures/'+document.getElementById('tree').value+'_provenance_'+document.getElementById('layout').value;document.getElementById('figure').src=stem+'.png';document.getElementById('svg').href=stem+'.svg';document.getElementById('png').href=stem+'.png';}</script></html>'''
    (P/'marker_provenance_review.html').write_text(doc,encoding='utf-8')
    print('COMPLETE: 14 plots, provenance table, PDF atlas and review page.',flush=True)

if __name__=='__main__':main()
