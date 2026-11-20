"""Update exact-cohort RF and ordination data without overwriting reviewed figures."""
import sys,json,inspect,statistics
from pathlib import Path
from worker import OUT,ROOT,load,save,sha
sys.path.insert(0,str(ROOT))
from terreus_tree_helpers import parse_newick,strong
DATA=ROOT/'comparative_visualization/2026-11-11/core66_combinations'
def main():
 accepted=OUT/'reports/core_support.validated.json'
 if not accepted.exists():
  save(OUT/'reports/new_core_comparison.json',{'status':'not_run','reason':'new core support did not pass'});return
 design=load(DATA/'design.json');keep=set(design['core_tips']);names=[c['name'] for c in design['combinations']]
 def splits(tree,core=False):
  if core:
   for t in tree.leaves():t.name='LOCAL_'+t.name.removesuffix('_transcripts').removesuffix('.transcripts')
  assert {t.name for t in tree.leaves()}==keep
  found=set()
  def visit(n):
   a={n.name} if not n.children else set().union(*(visit(c) for c in n.children));b=keep-a
   if n.children and len(a)>1 and len(b)>1:found.add(min((tuple(sorted(a)),tuple(sorted(b))),key=lambda x:(len(x),x)))
   return a
  visit(tree);return found
 def rf(a,b):return len(a^b)/(len(a)+len(b)) if a or b else 0.0
 newpath=OUT/'runs/core_support/core_GTR_F_R10.treefile'
 sets={'Core genome':splits(parse_newick(newpath.read_text()),True)}
 for name in names:sets[name]=splits(parse_newick((DATA/(name+'.treefile')).read_text()))
 labels=list(sets);matrix=[[rf(sets[a],sets[b]) for b in labels] for a in labels]
 rows=[]
 for name in names:
  values=[rf(splits(parse_newick(line)),sets['Core genome']) for line in (DATA/(name+'.ufboot')).read_text().splitlines() if line.strip()]
  assert len(values)==1000
  ordered=sorted(values)
  def quantile(q):
   x=q*(len(ordered)-1);a=int(x);b=min(a+1,len(ordered)-1);return ordered[a]+(x-a)*(ordered[b]-ordered[a])
  rows.append({'marker_choice':name,'RF':len(sets[name]^sets['Core genome']),'normalisation_denominator':len(sets[name])+len(sets['Core genome']),'normalised_RF':rf(sets[name],sets['Core genome']),'marker_bootstrap_replicates':1000,'marker_bootstrap_percentile_2_5':quantile(.025),'marker_bootstrap_percentile_97_5':quantile(.975),'interval_scope':'marker bootstrap versus fixed updated core ML topology'})
 save(OUT/'data/updated_core_RF.json',{'labels':labels,'matrix':matrix,'matched_tips':sorted(keep),'comparisons':rows,'core_tree_sha256':sha(newpath),'taxonomic_reference':'BenA + CaM + RPB2; topological similarity does not determine taxonomic authority'})
 # Report core bootstrap uncertainty separately, not as an unlabelled pooled interval.
 core_boot=[splits(parse_newick(line),True) for line in (OUT/'runs/core_support/core_GTR_F_R10.ufboot').read_text().splitlines() if line.strip()]
 save(OUT/'data/core_bootstrap_RF_to_marker_ML.json',{n:[rf(s,sets[n]) for s in core_boot] for n in names})
 try:
  import numpy as np,sklearn
  from sklearn.manifold import MDS
  params=inspect.signature(MDS).parameters
  kwargs={'n_components':2,'n_init':64,'max_iter':5000,'eps':1e-9,'random_state':20261118,'n_jobs':1}
  if 'metric_mds' in params:kwargs.update({'metric':'precomputed','metric_mds':False})
  else:kwargs.update({'metric':False,'dissimilarity':'precomputed'})
  if 'normalized_stress' in params:kwargs['normalized_stress']=True
  model=MDS(**kwargs);xy=model.fit_transform(np.asarray(matrix))
  save(OUT/'data/updated_core_NMDS.json',{'sklearn_version':sklearn.__version__,'arguments':kwargs,'stress':float(model.stress_),'stress_definition':'normalised stress' if 'normalized_stress' in kwargs else 'raw stress','rows':[{'tree':n,'x':float(xy[i,0]),'y':float(xy[i,1])} for i,n in enumerate(labels)]})
  nmds='PASS'
 except ImportError as e:nmds='DEPENDENCY_UNAVAILABLE: '+str(e)
 save(OUT/'reports/new_core_comparison.json',{'status':'RF_PASS','NMDS':nmds,'tips':len(keep),'combinations':len(names),'marker_bootstrap_RF_comparisons':7000,'core_bootstrap_comparisons':len(core_boot)*len(names),'figures':'reviewed previous figures preserved; new results are in this execution folder'})
if __name__=='__main__':main()
