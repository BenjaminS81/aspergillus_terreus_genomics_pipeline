from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from terreus_tree_helpers import parse_newick
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';D=P/'core66_combinations';F=P/'figures/core66_comparison'
keep=set(json.loads((D/'design.json').read_text())['core_tips']);ranking=json.loads((D/'ranking.json').read_text());genome=parse_newick(Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\0_4_mafft_iqtree\Sep28\concatenated.fa.treefile').read_text())
for t in genome.leaves():t.name='LOCAL_'+t.name.removesuffix('_transcripts').removesuffix('.transcripts')
def splits(root):
    ss=set()
    def visit(n):
        a={n.name} if not n.children else set().union(*(visit(c) for c in n.children));b=keep-a
        if n.children and len(a)>1 and len(b)>1:ss.add(min((tuple(sorted(a)),tuple(sorted(b))),key=lambda t:(len(t),t)))
        return a
    assert visit(root)==keep;return ss
target=splits(genome);allrows=[];stats=[];arrays=[]
for d in ranking:
    n=d['marker_choice'];values=[]
    with (D/(n+'.ufboot')).open() as f:
        for i,line in enumerate(f):
            if not line.strip():continue
            ss=splits(parse_newick(line));rf=len(ss^target)/(len(ss)+len(target));values.append(rf);allrows.append({'marker_choice':n,'bootstrap':i+1,'normalized_RF_to_core':rf})
    assert len(values)==1000;arrays.append(values);stats.append({'marker_choice':n,'replicates':len(values),'median_RF':float(np.median(values)),'percentile_2_5':float(np.quantile(values,.025)),'percentile_97_5':float(np.quantile(values,.975)),'ML_RF':d['normalized_RF_to_core']})
(D/'bootstrap_RF_summary.json').write_text(json.dumps(stats,indent=2))
with (D/'bootstrap_RF_values.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42});fig,ax=plt.subplots(figsize=(11,6));b=ax.boxplot(arrays,vert=False,patch_artist=True,showfliers=False,whis=(2.5,97.5),tick_labels=[x['marker_choice'].replace('_',' + ') for x in ranking])
for box in b['boxes']:box.set_facecolor('#c9e3df');box.set_edgecolor('#2a7d76')
ax.scatter([x['normalized_RF_to_core'] for x in ranking],range(1,8),color='#c17a45',s=35,zorder=5,label='Maximum-likelihood tree');ax.invert_yaxis();ax.set_xlim(0.45,1);ax.set_xlabel('Normalized RF to the fixed core-genome topology');ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,loc='lower right');ax.set_title('Marker topology uncertainty | 1,000 UFBoot trees per combination',loc='left',fontsize=16,weight='bold');fig.text(.06,.015,'Boxes: interquartile range; whiskers: 2.5th-97.5th percentiles. Core-tree uncertainty is not represented; replicate ranks are not paired across genes.',fontsize=9,color='#536778');fig.tight_layout(rect=[0,.05,1,1])
for ext in ['svg','pdf','png']:fig.savefig(F/('bootstrap_RF_distributions.'+ext),dpi=180,bbox_inches='tight')
plt.close(fig);print(json.dumps(stats,indent=2))
