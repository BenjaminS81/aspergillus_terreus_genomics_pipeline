from pathlib import Path
import json,csv,math,os,collections
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sklearn
from sklearn.manifold import MDS
from terreus_tree_helpers import parse_newick,strong
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis';F=P/'figures/tree_distances';F.mkdir(parents=True,exist_ok=True)
names=['BenA','CaM','RPB2','Concatenated'];trees={n:parse_newick((OLD/'trees'/((n if n!='Concatenated' else 'combined_all_strains')+'.treefile')).read_text()) for n in names}
core={'LOCAL_'+x['sample'] for x in json.loads((OLD/'genome_sources.json').read_text()) if x['in_core_tree']}
genome=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\0_4_mafft_iqtree\Sep28\concatenated.fa.treefile');gt=parse_newick(genome.read_text())
for t in gt.leaves():t.name='LOCAL_'+t.name.removesuffix('_transcripts').removesuffix('.transcripts')
assert {t.name for t in gt.leaves()}==core
def splitset(root,keep,mode):
    splits=set()
    def visit(n):
        desc=({n.name} & keep) if not n.children else set().union(*(visit(c) for c in n.children))
        if n.children and (mode=='all' or mode=='supported' and strong(n)):
            other=keep-desc
            if len(desc)>1 and len(other)>1:
                a,b=tuple(sorted(desc)),tuple(sorted(other));splits.add(min((a,b),key=lambda t:(len(t),t)))
        return desc
    visit(root);return splits
def incompatible(a,b,keep):
    a=set(a);b=set(b);return bool(a&b and a&(keep-b) and (keep-a)&b and (keep-a)&(keep-b))
# Independent four-tip checks include a polytomy and the root's duplicate split.
a=parse_newick('((A:1,B:1)100/100:1,(C:1,D:1)100/100:1);');b=parse_newick('((A:1,C:1)100/100:1,(B:1,D:1)100/100:1);');star=parse_newick('(A:1,B:1,C:1,D:1);');kk=set('ABCD')
assert len(splitset(a,kk,'all')^splitset(b,kk,'all'))==2 and len(splitset(a,kk,'all'))==1 and not splitset(star,kk,'all')
records=[];ordout=[];setsout={}
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans'})
cohorts=[('shared_reference_tips',set.intersection(*({t.name for t in root.leaves()} for root in trees.values())),trees),('core_66',core,trees),('core_66_with_genome',core,{**trees,'Core genome (exploratory)':gt})]
for cohort,keep,current in cohorts:
    nn=list(current);setsout[cohort]=sorted(keep);matrices={}
    for mode in ['all','supported']:
        if mode=='supported' and cohort.endswith('with_genome'):continue # genome run has no comparable support values
        ss={n:splitset(r,keep,mode) for n,r in current.items()};dist=np.zeros((len(nn),len(nn)))
        for i,n in enumerate(nn):
            for j,m in enumerate(nn):
                raw=len(ss[n]^ss[m]);den=len(ss[n])+len(ss[m]);dist[i,j]=raw/den if den else 0
                if i<j:records.append({'cohort':cohort,'shared_tips':len(keep),'mode':mode,'tree_1':n,'tree_2':m,'splits_1':len(ss[n]),'splits_2':len(ss[m]),'shared_splits':len(ss[n]&ss[m]),'RF':raw,'normalized_RF':dist[i,j],'splits_1_incompatible_with_tree_2':sum(any(incompatible(a,b,keep) for b in ss[m]) for a in ss[n]),'splits_2_incompatible_with_tree_1':sum(any(incompatible(a,b,keep) for b in ss[n]) for a in ss[m])})
        matrices[mode]=dist
        with (P/'analysis'/(cohort+'_'+mode+'_normalized_RF.tsv')).open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['tree']+nn);w.writerows([[n]+list(dist[i]) for i,n in enumerate(nn)])
        fig,ax=plt.subplots(figsize=(7.5,6));im=ax.imshow(dist,vmin=0,vmax=1,cmap='YlGnBu');ax.set_xticks(range(len(nn)),nn,rotation=30,ha='right');ax.set_yticks(range(len(nn)),nn)
        for i in range(len(nn)):
            for j in range(len(nn)):ax.text(j,i,f'{dist[i,j]:.3f}',ha='center',va='center',color='white' if dist[i,j]>.65 else '#163449')
        ax.set_title(f'Robinson-Foulds distance | {len(keep)} shared tips\n'+('All non-trivial splits' if mode=='all' else 'Only splits with SH-aLRT >=80 and UFBoot >=95'),loc='left',fontsize=12);fig.colorbar(im,ax=ax,label='RF / (splits in tree 1 + splits in tree 2)',shrink=.75);fig.tight_layout()
        for ext in ['svg','pdf','png']:fig.savefig(F/(cohort+'_'+mode+'_RF_heatmap.'+ext),dpi=170,bbox_inches='tight')
        plt.close(fig)
    # Ordination uses the full topology distance; weak-edge results remain a separate sensitivity panel.
    dist=matrices['all'];fig,axes=plt.subplots(1,2,figsize=(13,6))
    for metric,ax in zip([False,True],axes):
        model=MDS(n_components=2,metric_mds=metric,metric='precomputed',n_init=32,init='random',max_iter=3000,eps=1e-9,random_state=20261111,normalized_stress=True,n_jobs=1)
        xy=model.fit_transform(dist);method='Metric MDS' if metric else 'NMDS';err=float(model.stress_)
        for i,n in enumerate(nn):ax.scatter(*xy[i],s=100,color=['#287d8e','#cf8047','#7564a2','#2d947b','#757575'][i]);ax.annotate(n,xy[i],xytext=(6,7),textcoords='offset points',fontsize=10);ordout.append({'cohort':cohort,'method':method,'tree':n,'axis1':float(xy[i,0]),'axis2':float(xy[i,1]),'normalized_stress':err,'iterations':int(model.n_iter_)})
        ax.set_title(method+f' | Stress-1 = {err:.4f}',loc='left');ax.set_xlabel('Dimension 1 (arbitrary orientation)');ax.set_ylabel('Dimension 2');ax.spines[['top','right']].set_visible(False);ax.margins(.35);ax.set_aspect('equal',adjustable='datalim')
    fig.suptitle(f'Topology ordination | {len(keep)} shared tips',fontsize=16,weight='bold',x=.07,ha='left');fig.text(.07,.015,'Each point is a tree, not an isolate. Only four or five trees: this is descriptive, not a test of biological separation.',fontsize=10,color='#536778');fig.tight_layout(rect=[0,.06,1,.94])
    for ext in ['svg','pdf','png']:fig.savefig(F/(cohort+'_ordination.'+ext),dpi=170,bbox_inches='tight')
    plt.close(fig);print(cohort,'tips',len(keep),'RF',dist.round(3).tolist(),flush=True)
(P/'analysis/RF_pairwise_results.json').write_text(json.dumps(records,indent=2));(P/'analysis/ordination_results.json').write_text(json.dumps(ordout,indent=2));(P/'analysis/RF_shared_tip_sets.json').write_text(json.dumps(setsout,indent=2))
for filename,rows in [('RF_pairwise_results.tsv',records),('ordination_coordinates.tsv',ordout)]:
    with (P/'analysis'/filename).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(P/'analysis/tree_distance_methods.json').write_text(json.dumps({'RF':'Unrooted non-trivial bipartitions; symmetric difference after restriction to the cohort-wide shared tip set. Canonical split representation prevents double counting at the root. Normalization is RF/(observed splits in both trees).','support_filter':'Keep only original edges with both SH-aLRT >=80 and UFBoot >=95, then restrict splits to common tips. Polytomies retained. Resolution differences contribute to RF; incompatible-split counts are reported separately.','core_genome':'Exploratory original 4661-locus tree; no support values, so omitted from support-filtered comparison.','ordination':'scikit-learn '+sklearn.__version__+'; MDS(metric=precomputed, metric_mds=False/True), 2 dimensions, 32 random starts, 3000 iterations, eps=1e-9, seed=20261111, normalized_stress=True, one job.','limitations':'Four/five points; low stress alone does not validate topology or species boundaries. RF is insensitive to branch lengths and sensitive to arbitrary resolutions of near-identical sequences. Assembly aliases are retained as labelled in the supplied trees.'},indent=2))
