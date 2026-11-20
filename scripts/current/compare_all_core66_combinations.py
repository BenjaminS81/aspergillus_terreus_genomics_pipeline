from pathlib import Path
import os,json,csv,itertools
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from terreus_tree_helpers import parse_newick,strong,prune
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';D=P/'core66_combinations';F=P/'figures/core66_comparison';F.mkdir(parents=True,exist_ok=True)
design=json.loads((D/'design.json').read_text());keep=set(design['core_tips']);names=[x['name'] for x in design['combinations']];trees={n:parse_newick((D/(n+'.treefile')).read_text()) for n in names}
gpath=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\0_4_mafft_iqtree\Sep28\concatenated.fa.treefile');genome=parse_newick(gpath.read_text())
for t in genome.leaves():t.name='LOCAL_'+t.name.removesuffix('_transcripts').removesuffix('.transcripts')
trees={'Core genome':genome,**trees};assert all({t.name for t in r.leaves()}==keep for r in trees.values());nn=list(trees)
def splits(root,mode='all'):
    ss=set()
    def visit(n):
        s={n.name}&keep if not n.children else set().union(*(visit(c) for c in n.children));other=keep-s
        if n.children and len(s)>1 and len(other)>1 and (mode=='all' or mode=='supported' and strong(n) or mode=='length' and n.length>1.1e-6):ss.add(min((tuple(sorted(s)),tuple(sorted(other))),key=lambda x:(len(x),x)))
        return s
    visit(root);return ss
def rf(a,b):return len(a^b)/(len(a)+len(b)) if a or b else 0
def conflict(a,b):
    a=set(a);b=set(b);return bool(a&b and a&(keep-b) and (keep-a)&b and (keep-a)&(keep-b))
ss={n:splits(r) for n,r in trees.items()};lengthss={n:splits(r,'length') for n,r in trees.items()};ranks=[];pairs=[];matrix=np.zeros((8,8))
for i,n in enumerate(nn):
    for j,m in enumerate(nn):
        matrix[i,j]=rf(ss[n],ss[m])
        if i<j:pairs.append({'tree_1':n,'tree_2':m,'tips':66,'RF':len(ss[n]^ss[m]),'normalized_RF':matrix[i,j],'shared_splits':len(ss[n]&ss[m])})
for n in names:
    supported=splits(trees[n],'supported');shared=supported&ss['Core genome'];conflicting=sum(any(conflict(a,b) for b in ss['Core genome']) for a in supported)
    ranks.append({'marker_choice':n,'sites':next(d['sites'] for d in design['combinations'] if d['name']==n),'RF_to_core':len(ss[n]^ss['Core genome']),'normalized_RF_to_core':rf(ss[n],ss['Core genome']),'shared_core_splits':len(ss[n]&ss['Core genome']),'supported_marker_splits':len(supported),'supported_splits_matching_core':len(shared),'supported_splits_conflicting_core':conflicting,'supported_split_match_fraction':len(shared)/len(supported) if supported else None,'RF_after_short_edge_collapse':rf(lengthss[n],lengthss['Core genome'])})
ranks.sort(key=lambda d:(d['normalized_RF_to_core'],d['sites']))
for i,d in enumerate(ranks):d['rank_RF']=1+sum(x['normalized_RF_to_core']<d['normalized_RF_to_core'] for x in ranks)
(D/'ranking.json').write_text(json.dumps(ranks,indent=2));(D/'RF_pairwise.json').write_text(json.dumps(pairs,indent=2))
for fn,rows in [('ranking.tsv',ranks),('RF_pairwise.tsv',pairs)]:
    with (D/fn).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
with (D/'normalized_RF_matrix.tsv').open('w',newline='') as f:w=csv.writer(f,delimiter='\t');w.writerow(['tree']+nn);w.writerows([[n]+list(matrix[i]) for i,n in enumerate(nn)])
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans'})
def save(fig,name):
    for ext in ['svg','pdf','png']:fig.savefig(F/(name+'.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
fig,ax=plt.subplots(figsize=(11,9));im=ax.imshow(matrix,vmin=0,vmax=1,cmap='YlGnBu');labels=[x.replace('_',' + ') for x in nn];ax.set_xticks(range(8),labels,rotation=40,ha='right');ax.set_yticks(range(8),labels)
for i in range(8):
    for j in range(8):ax.text(j,i,f'{matrix[i,j]:.3f}',ha='center',va='center',color='white' if matrix[i,j]>.65 else '#213b54')
ax.set_title('All marker combinations versus the core genome\nExactly the same 66 tips; independently inferred marker trees',loc='left',fontsize=16,weight='bold');fig.colorbar(im,ax=ax,shrink=.7,label='Normalized unrooted Robinson-Foulds distance');fig.tight_layout();save(fig,'all_eight_RF_heatmap')
fig,axes=plt.subplots(1,2,figsize=(14,6),gridspec_kw={'width_ratios':[1.2,1]});ys=np.arange(7)
axes[0].barh(ys,[x['normalized_RF_to_core'] for x in ranks],color='#28877e');axes[0].set_yticks(ys,[x['marker_choice'].replace('_',' + ') for x in ranks]);axes[0].invert_yaxis();axes[0].set_xlim(0,1);axes[0].set_xlabel('Normalized RF to core genome (lower is closer)')
for i,d in enumerate(ranks):axes[0].text(d['normalized_RF_to_core']+.015,i,f"{d['normalized_RF_to_core']:.3f}",va='center')
axes[1].barh(ys,[x['supported_splits_matching_core'] for x in ranks],color='#28877e',label='Matches core split');axes[1].barh(ys,[x['supported_splits_conflicting_core'] for x in ranks],left=[x['supported_splits_matching_core'] for x in ranks],color='#c77851',label='Conflicts with core');axes[1].set_yticks(ys,[]);axes[1].invert_yaxis();axes[1].set_xlabel('Strong marker splits (SH-aLRT >=80; UFBoot >=95)');axes[1].legend(frameon=False)
for ax in axes:ax.spines[['top','right']].set_visible(False)
fig.suptitle('Which marker combination most closely reproduces the core topology?',fontsize=17,x=.06,ha='left',weight='bold');fig.text(.06,.015,'Topology similarity is not taxonomic validation. The original core tree has no branch-support analysis; short/weak marker splits can inflate RF.',fontsize=10,color='#536778');fig.tight_layout(rect=[0,.06,1,.92]);save(fig,'marker_combination_ranking')
fig,axes=plt.subplots(1,2,figsize=(15,7));ordinates=[];pal=['#253b55','#369297','#cc8650','#9275ae','#77a65a','#d05f79','#57a9c3','#c0aa41']
for metric,ax in zip([False,True],axes):
    model=MDS(n_components=2,metric='precomputed',metric_mds=metric,n_init=64,init='random',max_iter=5000,eps=1e-9,random_state=20261111,normalized_stress=True,n_jobs=1);xy=model.fit_transform(matrix);method='Metric MDS' if metric else 'NMDS'
    for i,n in enumerate(nn):
        ax.scatter(*xy[i],s=135 if i==0 else 80,color=pal[i],marker='*' if i==0 else 'o');ax.annotate(n.replace('_',' + '),xy[i],xytext=(7,7),textcoords='offset points',fontsize=9);ordinates.append({'method':method,'tree':n,'axis1':float(xy[i,0]),'axis2':float(xy[i,1]),'stress1':float(model.stress_),'iterations':int(model.n_iter_)})
    ax.set_title(method+f' | Stress-1 = {model.stress_:.4f}',loc='left',fontsize=14);ax.margins(.3);ax.set_xlabel('Dimension 1');ax.set_ylabel('Dimension 2');ax.spines[['top','right']].set_visible(False);ax.set_aspect('equal',adjustable='datalim')
fig.suptitle('Eight-tree topology space | exact 66-genome cohort',x=.06,ha='left',fontsize=19,weight='bold');fig.text(.06,.015,'Points are trees, not isolates. Distances are topological. Axis orientation is arbitrary; use the RF matrix for exact comparisons.',fontsize=10,color='#536778');fig.tight_layout(rect=[0,.055,1,.94]);save(fig,'all_eight_NMDS_MDS');(D/'ordination.json').write_text(json.dumps(ordinates,indent=2))
print(json.dumps(ranks,indent=2))
