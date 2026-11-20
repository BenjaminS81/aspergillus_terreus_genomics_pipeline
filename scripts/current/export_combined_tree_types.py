"""Two-page combined-tree export with prominent, provenance-aware type marks."""
import plot_marker_provenance as p
from pathlib import Path
import math,hashlib,json,csv
from matplotlib.lines import Line2D

OUT=p.P/'combined_type_marked';OUT.mkdir(exist_ok=True);(OUT/'qa').mkdir(exist_ok=True)
name='combined_all_strains';ann=p.annotations[name]
source=p.OLD/'sources/input_trees/combined_all_strains.treefile'
raw=source.read_bytes();root=p.parse_newick(raw.decode());tips=root.leaves()
clades=[x for x in p.allclades if x['tree']==name]
typed={k:a for k,a in ann.items() if a['type'] or a.get('genome_type_species')}
oldlabel=p.label
def label(k,a):
    return oldlabel(k,a).replace('[T]','[TYP]').replace('[TG]','[TYP-G]')
p.label=label
p.titles[name]='BenA + CaM + RPB2'
previous=p.P/'marker_trees_sequence_provenance.pdf';before=hashlib.sha256(previous.read_bytes()).hexdigest()
def save(fig,stem,pdf):
    circular=stem.endswith('_circular')
    if circular:
        ax=fig.axes[0]
        for i,t in enumerate(tips):
            if t.name in typed:ax.scatter(2*math.pi*i/len(tips),1.291,marker='*',s=70,facecolor='#e8af27',edgecolor='#533f12',linewidth=.45,zorder=9)
        text_axes=[ax]
    else:
        labels=fig.axes[1]
        for i,t in enumerate(tips):
            if t.name in typed:labels.scatter(-.022,i,marker='*',s=42,facecolor='#e8af27',edgecolor='#533f12',linewidth=.4,zorder=9,clip_on=False)
        text_axes=[labels]
    for ax in text_axes:
        for text in ax.texts:
            if '[TYP]' in text.get_text() or '[TYP-G]' in text.get_text():
                text.set_weight('bold');text.set_bbox({'boxstyle':'round,pad=0.12','facecolor':'#fff0bb','edgecolor':'none','alpha':.85})
    # Expand the existing type labels consistently in explanatory text.
    for text in fig.texts:
        s=text.get_text().replace('[T]','[TYP]').replace('[TG]','[TYP-G]')
        if s.startswith('BenA + CaM + RPB2 |'):s='BenA + CaM + RPB2 | Typstämme und Sequenzherkunft'
        text.set_text(s)
    fig.legend(handles=[Line2D([],[],marker='*',markerfacecolor='#e8af27',markeredgecolor='#533f12',lw=0,markersize=13,label='Typstamm: Goldstern + hervorgehobener Name')],loc='lower right',bbox_to_anchor=(.97,.025 if circular else .023),frameon=False,fontsize=11)
    fig.text(.035,.015 if circular else .041,'[TYP] 31 Referenz-Ex-Typ-Einträge   |   [TYP-G] 3 belegte Ex-Typ-Genome   |   Sequenzherkunft bleibt separat markiert.',fontsize=11,color='#533f12')
    pdf.savefig(fig,bbox_inches='tight',pad_inches=.3)
    suffix='kreisfoermig' if circular else 'linear'
    for ext in ('svg','png'):fig.savefig(OUT/('kombinierter_baum_typstaemme_'+suffix+'.'+ext),dpi=100,bbox_inches='tight',pad_inches=.3)
    p.plt.close(fig);p.gc.collect();print('EXPORTED',suffix,flush=True)
p.save=save
with p.PdfPages(OUT/'kombinierter_baum_typstaemme.pdf') as pdf:
    p.draw_circular(root,name,ann,clades,pdf)
    p.draw_linear(root,name,ann,pdf)
doc=p.pymupdf.open(OUT/'kombinierter_baum_typstaemme.pdf')
checks={'two_pages':len(doc)==2,'reference_type_entries':sum(bool(a['type']) for a in typed.values())==31,'genome_type_entries':sum(bool(a.get('genome_type_species')) for a in typed.values())==3,'previous_atlas_unchanged':hashlib.sha256(previous.read_bytes()).hexdigest()==before,'input_tree_unchanged':hashlib.sha256(source.read_bytes()).hexdigest()==hashlib.sha256(raw).hexdigest()}
for i,page in enumerate(doc):
    text=page.get_text();checks[f'page_{i+1}_type_labels']=text.count('[TYP]')==33 and text.count('[TYP-G]')==5
    page.get_pixmap(matrix=p.pymupdf.Matrix(.4,.4)).save(str(OUT/'qa'/f'page_{i+1}.png'))
    rect=page.rect
    if i==0:clip=p.pymupdf.Rect(rect.width*.55,rect.height*.42,rect.width*.99,rect.height*.8)
    else:clip=p.pymupdf.Rect(rect.width*.36,rect.height*.42,rect.width*.84,rect.height*.60)
    page.get_pixmap(matrix=p.pymupdf.Matrix(1,1),clip=clip).save(str(OUT/'qa'/f'page_{i+1}_detail.png'))
doc.close()
assert all(checks.values()),checks
(OUT/'validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
with (OUT/'markierte_typstaemme.tsv').open('w',newline='',encoding='utf-8') as h:
    writer=csv.DictWriter(h,fieldnames=['tip','taxon','strain','source','type_marker'],delimiter='\t');writer.writeheader()
    for k,a in typed.items():writer.writerow({'tip':k,'taxon':a['taxon'],'strain':a['strain'],'source':p.provenance(k),'type_marker':'TYP' if a['type'] else 'TYP-G'})
print('VALIDATED',checks,flush=True)
