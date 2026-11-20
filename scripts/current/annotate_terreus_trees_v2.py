from pathlib import Path
import json,csv,math,collections,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from terreus_tree_helpers import parse_newick,support,strong,prune
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis'
F=P/'figures/trees';F.mkdir(parents=True,exist_ok=True);OUT=P/'analysis';OUT.mkdir(exist_ok=True)
def js(p):return json.loads(p.read_text(encoding='utf-8'))
def short(s):return s.replace('Aspergillus_','A. ').replace('Aspergillus ','A. ').replace('_',' ')
meta={x['reference_id']:x for x in js(OLD/'tree_reference_metadata.json')};local=js(OLD/'local_metadata.json');tax={x['sample']:x for x in js(OLD/'sample_taxonomic_review.json')};sources=js(OLD/'genome_sources.json');occupancy={x['id']:x['available_loci'] for x in js(OLD/'tree_inputs/locus_occupancy.json')}
ncbi={x['summary']['sample_id']:x['summary'] for x in js(BASE/'review_reports/2026-11-10/ncbi_core_tree_metadata/core_tree_metadata_complete.json')}
annotations={}
for ident,d in meta.items():
    conflict=bool(d.get('accepted_type_name') and d['accepted_type_name']!=d['analysis_species']) or bool(d.get('historical_type_name'))
    annotations[ident]={'tip':ident,'taxon':short(d['analysis_species']),'strain':d['strain'],'role':'outgroup' if d['role']=='outgroup_reference' else 'type' if d['type_flag'] else 'historical' if ident.startswith('HIST2011_') else 'reference','species_ring':short(d['analysis_species']) if not conflict else 'Name conflict','taxonomy_evidence':'curated reference label' if not conflict else 'accepted/sequence-table or historical-name conflict','type':d['type_flag'],'type_name':short(d.get('historical_type_name') or d.get('accepted_type_name') or d['analysis_species']),'available_loci':occupancy.get(ident,0)}
for d in sources:
    s=d['sample'];ident='LOCAL_'+s;t=tax[s];old=local.get(ident,{}).get('original_labels',[])
    reported=''
    if old and '_' in old[0]:reported='A. '+old[0].split('_')[0]
    if not reported:
        n=ncbi.get(s,{}).get('ncbi_organism','').split(';')[0]
        reported=short(n) if n.startswith('Aspergillus ') else 'Aspergillus sp.'
    resolved=t['interpretation']=='three_locus_supported_concordance';name=short(t['concordant_supported_species']) if resolved else reported
    annotations[ident]={'tip':ident,'taxon':name,'strain':s,'role':'core' if d['in_core_tree'] else 'additional','species_ring':name if resolved else 'Unresolved','taxonomy_evidence':'three-locus supported type placement' if resolved else 'reported label; species unresolved','type':False,'type_name':'','available_loci':occupancy.get(ident,0)}
genome_types={x['sample']:x for x in js(OUT/'genome_type_anchor_audit.json')['accepted_genome_type_anchors']}
for k,d in annotations.items():
    anchor=genome_types.get(k.removeprefix('LOCAL_')) if k.startswith('LOCAL_') else None
    d['genome_type_species']=short(anchor['species']) if anchor else ''
    d['genome_type_culture']=anchor['culture'] if anchor else ''
    if anchor:d['taxon']=short(anchor['species']);d['species_ring']=d['taxon'];d['taxonomy_evidence']='documented ex-type genome culture'
all_species=sorted({d['species_ring'] for d in annotations.values()}-{'Unresolved','Name conflict'})
cmap=plt.get_cmap('turbo');colors={s:matplotlib.colors.to_hex(cmap(.05+.9*i/max(1,len(all_species)-1))) for i,s in enumerate(all_species)};colors.update({'Unresolved':'#c9cdd2','Name conflict':'#a7a1aa','Mixed type species':'#b6b7bf','No supported type clade':'#f3f4f6'})
rolecolors={'core':'#147f78','additional':'#4c9fbc','type':'#203b61','reference':'#8194aa','historical':'#b0b8c1','outgroup':'#c28c2e'}
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'font.size':9})
def label(k):
    d=annotations[k];q=' *' if d['taxonomy_evidence']=='reported label; species unresolved' else ''
    t=' [T]' if d['type'] else ''
    if d['species_ring']=='Name conflict':t+=' [name flag]'
    if d.get('genome_type_species'):t+=' [TG]'
    return d['taxon']+q+' | '+d['strain']+t+(' [2011]' if k.startswith('HIST2011_') else '')
def index_tree(root):
    leaves={};paths={}
    def visit(n,path):
        if not n.children:leaves[id(n)]=[n.name];paths[n.name]=path+[n]
        else:
            for c in n.children:visit(c,path+[n])
            leaves[id(n)]=sum((leaves[id(c)] for c in n.children),[])
    visit(root,[]);return leaves,paths
def type_clades(root,name):
    leaves,paths=index_tree(root);types={k:d['analysis_species'] for k,d in meta.items() if d['type_flag']};type_sets={i:{types[k] for k in names if k in types} for i,names in leaves.items()}
    chosen={};rows=[];selected={}
    for k,path in paths.items():
        pure=[n for n in reversed(path[:-1]) if len(type_sets[id(n)])==1 and strong(n)]
        mixed=[n for n in reversed(path[:-1]) if 1<len(type_sets[id(n)])<=8 and strong(n)]
        n=pure[-1] if pure else mixed[0] if mixed else None
        if n is not None:selected[id(n)]=n
        chosen[k]=n
    ordered=sorted(selected.values(),key=lambda n:(sorted(type_sets[id(n)]),len(leaves[id(n)])))
    ids={id(n):'C'+str(i+1).zfill(2) for i,n in enumerate(ordered)}
    clades=[]
    for n in ordered:
        tt=sorted(type_sets[id(n)]);a,b=support(n);clades.append({'tree':name,'clade_id':ids[id(n)],'type_species':'; '.join(short(x) for x in tt),'type_reference_tips':';'.join(k for k in leaves[id(n)] if k in types),'SH_aLRT':a,'UFBoot':b,'descendant_count':len(leaves[id(n)]),'descendant_tips':leaves[id(n)],'status':'single_type_species' if len(tt)==1 else 'mixed_type_species'})
    for k,n in chosen.items():
        tt=sorted(type_sets[id(n)]) if n is not None else [];a,b=support(n) if n is not None else (None,None)
        rows.append({'tree':name,**annotations[k],'display_label':label(k),'type_clade':ids[id(n)] if n is not None else 'unresolved','type_clade_species':'; '.join(short(x) for x in tt),'type_clade_category':short(tt[0]) if len(tt)==1 else 'Mixed type species' if tt else 'No supported type clade','SH_aLRT':a,'UFBoot':b,'genome_type_clade_category':annotations[k].get('genome_type_species') or 'No supported type clade'})
    return {x['tip']:x for x in rows},clades
def positions(root,clado=False):
    tips=root.leaves();ys={id(n):i for i,n in enumerate(tips)};xs={id(root):0}
    def visit(n):
        for c in n.children:xs[id(c)]=xs[id(n)]+(1 if clado else c.length);visit(c)
        if n.children:ys[id(n)]=sum(ys[id(c)] for c in n.children)/len(n.children)
    visit(root);return xs,ys,tips
def save(fig,stem):
    for ext in ['svg','pdf','png']:fig.savefig(F/(stem+'.'+ext),dpi=145,bbox_inches='tight',pad_inches=.3)
    plt.close(fig)
def linear(root,title,stem,ann,clades,projected=False):
    xs,ys,tips=positions(root);mx=max(xs.values());n=len(tips);fig=plt.figure(figsize=(24,max(10,n*.18+3)))
    ax=fig.add_axes([.035,.045,.38,.9]);lab=fig.add_axes([.43,.045,.39,.9],sharey=ax);tracks=fig.add_axes([.83,.045,.15,.9],sharey=ax)
    def draw(node):
        if node.children:
            ax.plot([xs[id(node)]]*2,[min(ys[id(c)] for c in node.children),max(ys[id(c)] for c in node.children)],c='#687888',lw=.65)
            for c in node.children:ax.plot([xs[id(node)],xs[id(c)]],[ys[id(c)]]*2,c='#687888',lw=.65);draw(c)
            if strong(node):ax.scatter(xs[id(node)],ys[id(node)],s=4,color='#243c55')
        else:ax.plot([xs[id(node)],mx*1.01],[ys[id(node)]]*2,c='#d7dde3',lw=.4,ls=':')
    draw(root)
    for i,t in enumerate(tips):
        d=ann[t.name];lab.text(0,i,label(t.name),va='center',fontsize=8,color=rolecolors[d['role']],weight='bold' if d['role'] in ['core','additional'] else 'normal')
        for j,col in enumerate([rolecolors[d['role']],{1:'#d7794c',2:'#e6ad51',3:'#408f87'}.get(d['available_loci'],'#eee'),colors[d['species_ring']],colors[d['type_clade_category']],colors.get(d.get('genome_type_clade_category','No supported type clade'),'#f3f4f6')]):tracks.add_patch(plt.Rectangle((j-.38,i-.47),.76,.94,color=col,lw=0))
        tracks.text(4.65,i,d['type_clade'],va='center',fontsize=7,color='#30455a')
    for a in [ax,lab,tracks]:a.set_ylim(n+.5,-2)
    ax.set_xlim(0,mx*1.02);ax.set_yticks([]);ax.spines[['top','right','left']].set_visible(False);ax.set_xlabel('Substitutions per site (uncompressed)')
    lab.set_xlim(0,1);lab.axis('off');tracks.set_xlim(-.7,6.8);tracks.set_xticks(range(5),['Role','Loci','Species','Ref. type clade','Genome type'],rotation=45,ha='right');tracks.xaxis.tick_top();tracks.set_yticks([]);tracks.spines[:].set_visible(False)
    fig.suptitle(title,fontsize=22,x=.035,y=.985,ha='left',weight='bold',color='#213b54')
    fig.text(.035,.965,f'{n} tips | taxon + strain labels | dots: SH-aLRT >=80 and UFBoot >=95 | * reported species, unresolved by the three markers',fontsize=10,color='#526778')
    fig.text(.035,.015,('Reference ring is projected; outer genome-type ring marks [TG] cultures and compatible supported clades within this tree.' if projected else 'Type-clade codes are defined in the key. [TG] marks a documented genome ex-type culture; the outer ring locates these genomes.'),fontsize=10,color='#526778');save(fig,stem+'_linear')
def circular(root,title,stem,ann,clades,projected=False):
    xs,ys,tips=positions(root,True);n=len(tips);mx=max(xs.values());theta={id(t):2*math.pi*i/n for i,t in enumerate(tips)}
    def calc(node):
        for c in node.children:calc(c)
        if node.children:theta[id(node)]=sum(theta[id(c)] for c in node.children)/len(node.children)
    calc(root);size=34 if n>160 else 23;fig,ax=plt.subplots(figsize=(size,size),subplot_kw={'projection':'polar'});step=2*math.pi/n
    # Pale sectors mark maximal supported single-type clades, not provisional species calls.
    single=sorted([c for c in clades if c['status']=='single_type_species'],key=lambda x:-len(x['descendant_tips']));painted=set()
    idx={t.name:i for i,t in enumerate(tips)}
    for c in single:
        inds=[idx[k] for k in c['descendant_tips'] if k in idx]
        if not inds or len(inds)<2 or set(inds)&painted:continue
        lo,hi=min(inds),max(inds)
        if hi-lo+1!=len(inds):continue
        ax.bar((lo+hi)*step/2,1.02,width=(hi-lo+1)*step,bottom=0,color=colors[c['type_species']],alpha=.055,edgecolor='none',zorder=0);painted.update(inds)
    def draw(node):
        r=xs[id(node)]/mx
        if node.children:
            lo=min(theta[id(c)] for c in node.children);hi=max(theta[id(c)] for c in node.children);ax.plot([lo+(hi-lo)*i/50 for i in range(51)],[r]*51,c='#70808f',lw=.65)
            for c in node.children:ax.plot([theta[id(c)]]*2,[r,xs[id(c)]/mx],c='#70808f',lw=.65);draw(c)
            if strong(node):ax.scatter(theta[id(node)],r,s=4,c='#253e59',zorder=3)
        else:
            a=theta[id(node)];d=ann[node.name];ax.plot([a,a],[r,1.03],c='#c9d3db',lw=.45,ls=':')
            cols=[rolecolors[d['role']],{1:'#d7794c',2:'#e6ad51',3:'#408f87'}.get(d['available_loci'],'#eee'),colors[d['species_ring']],colors[d['type_clade_category']],colors.get(d.get('genome_type_clade_category','No supported type clade'),'#f3f4f6')]
            for j,col in enumerate(cols):ax.bar(a,.034,width=step*.94,bottom=1.045+j*.045,color=col,lw=0,zorder=3)
            deg=math.degrees(a);left=90<deg<270;txt=label(node.name)
            ax.text(a,1.31,txt,rotation=deg+(180 if left else 0),rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=8 if n>160 else 9,color=rolecolors[d['role']],weight='bold' if d['role'] in ['core','additional'] else 'normal')
    draw(root)
    # Clade codes sit just inside the outer rings; exact taxa and supports are in the key.
    runs=[]
    for i,t in enumerate(tips):
        code=ann[t.name]['type_clade']
        if runs and runs[-1][2]==code:runs[-1][1]=i
        else:runs.append([i,i,code])
    for lo,hi,code in runs:
        if code=='unresolved' or hi-lo<2:continue
        a=(lo+hi)*step/2;deg=math.degrees(a);ax.text(a,1.01,code,fontsize=8,color='#213b54',ha='center',va='center',rotation=deg-90+(180 if 90<deg<270 else 0),rotation_mode='anchor')
    ax.set_ylim(0,1.96);ax.axis('off');fig.suptitle(title,x=.035,y=.98,ha='left',fontsize=27,weight='bold',color='#213b54')
    fig.text(.035,.95,f'{n} tips | circular cladogram: branch lengths not represented; uncompressed linear companion supplied\nTaxon + strain at every leaf. * Reported species, unresolved by the three markers. [T] Ex-type reference.',fontsize=12,color='#526778')
    fig.text(.035,.035,('Rings: role | available loci | reviewed species | projected reference type clade | genome-type evidence\n[TG] = documented genome ex-type culture. Outer ring: these anchors and compatible supported clades in this tree.\nReference placement is projected; the three genome type cultures provide limited in-tree taxon coverage.' if projected else 'Rings: role | available loci | reviewed species | supported reference type clade | genome-type evidence\n[T] = reference ex-type; [TG] = documented genome ex-type culture. Grey species: unresolved; lavender: name conflict.\nGrey type ring: mixed types; white: no supported type clade. Clade codes refer to the tree-specific key.'),fontsize=12,color='#526778')
    fig.legend(handles=[Patch(color=c,label=s) for s,c in rolecolors.items()],loc='lower right',bbox_to_anchor=(.97,.025),ncol=3,frameon=False,fontsize=11)
    fig.subplots_adjust(left=.02,right=.98,top=.92,bottom=.075);save(fig,stem+'_circular')

def main():
    allrows=[];allclades=[];trees={}
    for name,title in [('BenA','BenA'),('CaM','Calmodulin (CaM)'),('RPB2','RNA polymerase II (RPB2)'),('combined_all_strains','Concatenated BenA + CaM + RPB2')]:
        root=parse_newick((OLD/'trees'/(name+'.treefile')).read_text());ann,clades=type_clades(root,name);allrows+=list(ann.values());allclades+=clades
        trees[name]={'tips':len(root.leaves()),'source':str(OLD/'trees'/(name+'.treefile')),'sha256':hashlib.sha256((OLD/'trees'/(name+'.treefile')).read_bytes()).hexdigest()}
        for focused in [False,True]:
            keep={k for k,d in annotations.items() if k.startswith('LOCAL_') or d['type'] or (k in meta and (meta[k].get('accepted_type_name') or meta[k].get('historical_type_name')))}
            r=prune(root,keep) if focused else root;stem=name+('_type_focused' if focused else '_all_references');tt=title+(' | type-focused display' if focused else ' | all reference records')
            linear(r,tt,stem,ann,clades);circular(r,tt,stem,ann,clades)
        print('ANNOTATED',name,len(root.leaves()),'tips',len(clades),'supported type-clade entries',flush=True)
    (OUT/'tree_annotations.json').write_text(json.dumps(allrows,indent=2));(OUT/'type_clades.json').write_text(json.dumps(allclades,indent=2));(OUT/'taxon_palette.json').write_text(json.dumps(colors,indent=2));(OUT/'sample_annotations.json').write_text(json.dumps(annotations,indent=2));(OUT/'tree_sources.json').write_text(json.dumps(trees,indent=2))
    for filename,rows in [('tree_annotations.tsv',allrows),('type_clades.tsv',[{**x,'descendant_tips':';'.join(x['descendant_tips'])} for x in allclades])]:
        with (OUT/filename).open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(figsize=(13,10));ax.axis('off');ax.set_title('Shared species colours and annotation policy',loc='left',fontsize=19,weight='bold',color='#213b54')
    handles=[Patch(color=c,label=s) for s,c in colors.items()];ax.legend(handles=handles,loc='upper left',ncol=2,frameon=False,fontsize=10)
    fig.text(.07,.075,'Species ring: curated reference species, or project species supported independently by all three loci.\nNo project A. terreus sensu-stricto assignment is forced where these markers leave its boundary unresolved.\nNo formal subspecies annotations were established in the source metadata.\nType-clade ring: tree-specific ex-type anchors; both SH-aLRT >=80 and UFBoot >=95 required.\nType-focused figures prune the full trees for display; they are not independently inferred trees.',fontsize=10,color='#526778');save(fig,'annotation_legend')

if __name__=="__main__":main()
