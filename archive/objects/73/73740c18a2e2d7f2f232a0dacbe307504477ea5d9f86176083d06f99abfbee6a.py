"""Read-only source audit; new, bounded clade-congruence release (no tree inference)."""
from pathlib import Path
import os, sys, json, hashlib, itertools, collections, csv, math, copy, shutil
ROOT=Path(__file__).resolve().parent
from circos_review_common import constrain_process, fasta
LIMITS=constrain_process()
from terreus_tree_helpers import Node, parse_newick, prune, strong
import numpy as np
OUT=ROOT/'stakeholder_release/2026-11-20_v4_clade_congruence_cgmlst'
for d in ('data','figures','qa','scripts'):(OUT/d).mkdir(parents=True,exist_ok=True)
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(name,v):(OUT/'data'/name).write_text(json.dumps(v,indent=2),encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ANN_PATH=ROOT/'stakeholder_release/2026-11-20_v3_original_taxon_labels/data/core_tip_assignments.json'
ANN={a['tip']:a for a in read(ANN_PATH)};NAMES=sorted(ANN);ALL=set(NAMES)
P=ROOT/'comparative_visualization/2026-11-11/core66_combinations'
NEW=ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene'
CORE=ROOT/'stakeholder_release/2026-11-20_v1_type_monophyly/data/core66_display_rooted.treefile'
KEYS=['BenA','CaM','RPB2','BenA_CaM','BenA_RPB2','CaM_RPB2','BenA_CaM_RPB2','Three_gene_reference230']
PATHS={'Core':CORE,**{k:P/(k+'.treefile') for k in KEYS[:-1]},KEYS[-1]:NEW/'runs/complete230.treefile'}
INPUTS={str(p):sha(p) for p in [*PATHS.values(),ANN_PATH,P/'design.json',P/'execution.json',NEW/'status.json']}
def nodes(r):
 yield r
 for c in r.children:yield from nodes(c)
def canon(s,univ):
 s=frozenset(s);o=frozenset(univ-s)
 return min((s,o),key=lambda a:(len(a),tuple(sorted(a))))
def splitmap(r,keep=None):
 univ=set(t.name for t in r.leaves()) if keep is None else set(keep);result={}
 def rec(n):
  s=set().union(*(rec(c) for c in n.children)) if n.children else ({n.name}&univ)
  if n is not r and 1<len(s)<len(univ)-1:
   result.setdefault(canon(s,univ),[]).append(n)
  return s
 assert rec(r)==univ
 return result
def rf(a,b):
 sa=set(splitmap(a));sb=set(splitmap(b));return len(sa^sb)/(len(sa)+len(sb)) if sa or sb else None
def distances(r,unit=False):
 paths={}
 def rec(n,path,dep):
  path=path+[(id(n),dep)]
  if not n.children:paths[n.name]=path
  for c in n.children:rec(c,path,dep+(1 if unit else c.length))
 rec(r,[],0);names=sorted(paths);mat=np.zeros((len(names),len(names)))
 for i,a in enumerate(names):
  da=dict(paths[a])
  for j in range(i+1,len(names)):
   b=names[j];lca=next(d for n,d in reversed(paths[b]) if n in da)
   mat[i,j]=mat[j,i]=paths[a][-1][1]+paths[b][-1][1]-2*lca
 return mat
def reroot(r,tip='LOCAL_SRR10507650'):
 ns=list(nodes(r));lookup={id(n):n for n in ns};adj={id(n):[] for n in ns}
 for n in ns:
  for c in n.children:
   lab=c.name if c.children else ''
   adj[id(n)].append((id(c),c.length,lab));adj[id(c)].append((id(n),c.length,lab))
 u=id(next(t for t in r.leaves() if t.name==tip));v,length,label=adj[u][0]
 def walk(x,prev,l,lab):
  cs=[walk(y,x,z,w) for y,z,w in adj[x] if y!=prev]
  n=Node(lookup[x].name if not cs else lab,l,cs)
  if len(cs)==1:cs[0].length+=l;return cs[0]
  return n
 return Node('',0,[walk(u,v,length/2,label),walk(v,u,length/2,label)])
def newick(n):return ('('+','.join(newick(c) for c in n.children)+')' if n.children else '')+n.name+':'+format(n.length,'.14g')
def layout_order(n,rank=None):
 for c in n.children:layout_order(c,rank)
 if rank is None:n.children.sort(key=lambda c:(-len(c.leaves()),min(t.name for t in c.leaves())))
 elif len(n.children)>1:
  assert len(n.children)<=6,'Review high-degree display polytomy'
  seq={id(c):[rank[t.name] for t in c.leaves()] for c in n.children}
  def cost(order):return sum(a>b for i,c in enumerate(order) for d in order[i+1:] for a in seq[id(c)] for b in seq[id(d)])
  n.children=list(min(itertools.permutations(n.children),key=cost))
def quartet_codes(d,q):
 a,b,c,e=q.T;s=np.stack((d[a,b]+d[c,e],d[a,c]+d[b,e],d[a,e]+d[b,c]),axis=1)
 m=s.min(axis=1);codes=s.argmin(axis=1).astype(np.int8);codes[(np.isclose(s,m[:,None],rtol=0,atol=1e-8)).sum(axis=1)!=1]=-1
 return codes
# Biological algorithm checks: identical quartet, changed split, unresolved star, pruning and reroot invariance.
t1=parse_newick('((a:1,b:1):1,(c:1,d:1):1);');t2=parse_newick('((a:1,c:1):1,(b:1,d:1):1);');star=parse_newick('(a:1,b:1,c:1,d:1);');q=np.array([[0,1,2,3]])
assert rf(t1,t1)==0 and rf(t1,t2)==1
assert quartet_codes(distances(t1,True),q)[0]==0 and quartet_codes(distances(t2,True),q)[0]==1
assert quartet_codes(distances(star,True),q)[0]==-1
assert np.allclose(distances(t1),distances(reroot(t1,'a')))
TREES={};VALID=[]
for k,p in PATHS.items():
 raw=parse_newick(p.read_text());tips=[t.name for t in raw.leaves()];assert len(tips)==len(set(tips))
 assert ALL<=set(tips)
 r=prune(raw,ALL);r.length=0;before=set(splitmap(r));d=distances(r)
 r=reroot(r);assert before==set(splitmap(r));delta=float(np.max(np.abs(d-distances(r))));assert delta<1e-9
 TREES[k]=r;VALID.append({'tree':k,'source_tips':len(tips),'matched_tips':len(r.leaves()),'unique_matching_ids':set(t.name for t in r.leaves())==ALL,'preserved_unrooted_splits':True,'max_patristic_change':delta})
layout_order(TREES['Core']);ORDER=[t.name for t in TREES['Core'].leaves()];rank={t:i for i,t in enumerate(ORDER)}
for k in KEYS:layout_order(TREES[k],rank)
for k,r in TREES.items():(OUT/'data'/(k+'.matched66.treefile')).write_text(newick(r)+';')
oldseq=dict(fasta(P/'BenA_CaM_RPB2.fasta'));newseq=dict(fasta(NEW/'inputs/complete230.fasta'))
seqdiff=[s for s in NAMES if oldseq[s]!=newseq[s]]
assert not seqdiff,'Old marker alignments differ from current three-gene source; review versions'
save('input_sequence_equivalence.json',{'identical_aligned_sequences':66,'sites':2076,'partitions':{'BenA':[1,515],'CaM':[516,1062],'RPB2':[1063,2076]},'different_tips':seqdiff,'scope':'Seven 66-tip inference inputs versus matching rows of complete230; tree inference and taxon sampling differ.'})
for p in [P/'BenA_CaM_RPB2.fasta',NEW/'inputs/complete230.fasta']:INPUTS[str(p)]=sha(p)
GROUPS={s:{t for t in ALL if ANN[t]['assigned_species']==s} for s in sorted({a['assigned_species'] for a in ANN.values()})}
GROUPS['terreus + montoensis']=GROUPS['A. terreus']|GROUPS['A. montoensis']
BOOT={};BTCOUNTS={}
for k in KEYS:
 p=(P/(k+'.ufboot')) if k!=KEYS[-1] else NEW/'runs/complete230.ufboot';INPUTS[str(p)]=sha(p)
 count=collections.Counter();n=0
 with p.open() as f:
  for line in f:
   if not line.strip():continue
   bt=parse_newick(line.strip());assert ALL<=set(t.name for t in bt.leaves());count.update(splitmap(bt,ALL).keys());n+=1
 assert n==1000,(k,n);BOOT[k]=count;BTCOUNTS[k]=n;print('BOOTSTRAPS',k,n,flush=True)
SPLITS={k:splitmap(r) for k,r in TREES.items()};coreS=set(SPLITS['Core']);clades=[]
for g,target in GROUPS.items():
 for k,r in TREES.items():
  n=len(target);key=canon(target,ALL);mono=None if n<2 else key in SPLITS[k]
  candidates=[s for sp in SPLITS[k] for s in (set(sp),ALL-set(sp))]
  def score(s):return (2*len(s&target)/(len(s)+n),len(s&target)/len(s),-len(s))
  best=max(candidates,key=score) if n>1 else set(target);tp=best&target
  subRF=rf(prune(TREES['Core'],target),prune(r,target)) if n>=4 else None
  clades.append({'group':g,'n':n,'tree':k,'strict_split_recovered':mono,'bootstrap_split_frequency_pct':100*BOOT[k][key]/BTCOUNTS[k] if n>=2 and k in BOOT else None,'within_group_normalized_RF_to_core':subRF,'nearest_split_F1':2*len(tp)/(len(best)+n) if n>=2 else None,'nearest_split_precision':len(tp)/len(best) if n>=2 else None,'nearest_split_recall':len(tp)/n if n>=2 else None,'nearest_split_extra_tips':sorted(best-target),'nearest_split_missing_tips':sorted(target-best)})
save('clade_congruence.json',clades)
def incompatible(a,b):return bool(a&b and a&(ALL-b) and (ALL-a)&b and (ALL-a)&(ALL-b))
edges=[]
for i,s in enumerate(sorted(coreS,key=lambda x:(len(x),tuple(sorted(x)))),1):
 e={'edge_id':f'CG{i:02d}','smaller_side':sorted(s),'trees':{}}
 for k in KEYS:
  conflict=[sp for sp,ns in SPLITS[k].items() if incompatible(set(s),set(sp)) and k!=KEYS[-1] and any(strong(n) for n in ns)]
  e['trees'][k]={'exact_split_present':s in SPLITS[k],'bootstrap_split_frequency_pct':100*BOOT[k][s]/1000,'strong_marker_conflicting_splits':len(conflict) if k!=KEYS[-1] else None}
 edges.append(e)
save('core_edge_congruence.json',edges)
# Exhaustive, topological quartet comparison. Streaming batches avoid large temporary arrays.
D={k:distances(r,True) for k,r in TREES.items()};QT={k:collections.Counter() for k in KEYS};TIPNUM={k:np.zeros(66,dtype=int) for k in KEYS};TIPDEN={k:np.zeros(66,dtype=int) for k in KEYS}
iterator=itertools.combinations(range(66),4);total=0
while True:
 batch=list(itertools.islice(iterator,20000))
 if not batch:break
 q=np.asarray(batch,dtype=np.int16);total+=len(q);c=quartet_codes(D['Core'],q)
 for k in KEYS:
  z=quartet_codes(D[k],q);both=(c>=0)&(z>=0);same=both&(c==z);diff=both&(c!=z)
  QT[k].update({'both_resolved':int(both.sum()),'agree':int(same.sum()),'disagree':int(diff.sum()),'only_core_resolved':int(((c>=0)&(z<0)).sum()),'only_marker_resolved':int(((c<0)&(z>=0)).sum()),'both_unresolved':int(((c<0)&(z<0)).sum())})
  TIPNUM[k]+=np.bincount(q[diff].ravel(),minlength=66);TIPDEN[k]+=np.bincount(q[both].ravel(),minlength=66)
assert total==math.comb(66,4)==720720
summary=[]
for k in KEYS:
 item={'tree':k,'normalized_RF':rf(TREES['Core'],TREES[k]),'shared_core_splits':len(coreS&set(SPLITS[k])),'total_quartets':total,**QT[k],'quartet_agreement_pct':100*QT[k]['agree']/QT[k]['both_resolved']}
 # Near-zero-edge sensitivity, separate from the literal topology results.
 def positive(n):return any(x.length>1.1e-6 for x in n)
 sc={s for s,ns in SPLITS['Core'].items() if positive(ns)};sm={s for s,ns in SPLITS[k].items() if positive(ns)}
 item['normalized_RF_edges_over_1_1e_6']=len(sc^sm)/(len(sc)+len(sm));summary.append(item)
save('global_comparisons.json',summary)
save('tip_quartet_discordance.json',[{'tip':t,'species':ANN[t]['assigned_species'],'tree':k,'discordant':int(TIPNUM[k][i]),'both_resolved':int(TIPDEN[k][i]),'discordance_pct':float(100*TIPNUM[k][i]/TIPDEN[k][i])} for i,t in enumerate(NAMES) for k in KEYS])
# Existing protein/transcript-model copy counts: candidate inventory only, no allele calling.
MATRIX=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\0_3_orthofinder\terrei3\Results_Sep28\Orthogroups\Orthogroups.GeneCount.tsv')
rows=[];ids=[]
with MATRIX.open(newline='') as f:
 reader=csv.DictReader(f,delimiter='\t');cols=[c for c in reader.fieldnames if c not in ('Orthogroup','Total')]
 assert {'LOCAL_'+c.removesuffix('_transcripts') for c in cols}==ALL
 for row in reader:ids.append(row['Orthogroup']);rows.append([int(row[c]) for c in cols])
arr=np.array(rows,dtype=np.int16);colnames=['LOCAL_'+c.removesuffix('_transcripts') for c in cols];INPUTS[str(MATRIX)]=sha(MATRIX)
canonrows=read(ROOT/'gate_preparation/2026-11-18_v1/data/canonical_assemblies.json');dedup={s:r['sequence_sha256'] for r in canonrows for s in r['samples']}
inventory=[];candidate_lists={}
panels={'Current 66-genome panel':ALL,**{g:s for g,s in GROUPS.items() if '+' not in g}}
for g,tips in panels.items():
 ix=[colnames.index(t) for t in sorted(tips)];a=arr[:,ix];zero_dups=(a<=1).all(axis=1);strict=(a==1).all(axis=1);s95=zero_dups&((a==1).sum(axis=1)>=math.ceil(.95*len(ix)))
 unique={dedup.get(t.removeprefix('LOCAL_'),'UNRESOLVED:'+t) for t in tips}
 inv={'panel':g,'tip_labels':len(ix),'distinct_sequence_hashes_or_unresolved_ids':len(unique),'all_single_copy_100pct':int(strict.sum()),'at_least_95pct_single_copy_no_duplicates':int(s95.sum()),'orthogroups_with_any_multiple_copy_model':int((a>1).any(axis=1).sum()),'status':'candidate inventory only; annotation, paralogy, CDS integrity, provenance and external validation gates remain'}
 inventory.append(inv);candidate_lists[g]={'strict_100pct':[ids[i] for i in np.where(strict)[0]],'occupancy95pct_no_duplicates':[ids[i] for i in np.where(s95)[0]]}
corelist=MATRIX.parent/'core.txt';origcore=set(corelist.read_text().split());strictset=set(candidate_lists['Current 66-genome panel']['strict_100pct'])
save('cgmlst_candidate_inventory.json',{'matrix_path':str(MATRIX),'orthogroups':len(ids),'panels':inventory,'existing_core_list_count':len(origcore),'existing_core_vs_count_strict':{'intersection':len(origcore&strictset),'old_core_not_strict':sorted(origcore-strictset),'strict_not_old_core':sorted(strictset-origcore)},'caveat':'These are copy-count screens over current gene models, not validated cgMLST loci. This panel represents nine provisionally assigned species, not all section Terrei diversity. Singleton and n=2 panels cannot establish species-wide core status.'})
save('cgmlst_candidate_orthogroups.json',candidate_lists)
save('core_tip_assignments.json',list(ANN.values()))
save('analysis_manifest.json',{'resource_limits':LIMITS,'tree_keys':list(TREES),'comparisons':8,'matched_tips':66,'bootstrap_replicates':BTCOUNTS,'all_quartets_per_comparison':total,'core_order':ORDER,'input_hashes':INPUTS,'validations':VALID,'method':'Unrooted canonical splits; branch-length-preserving pruning and display rerooting; exact same 66 IDs; all 720720 quartets from unit-edge tree distances; 1000 UFBoot trees per marker; deterministic rotations minimize crossings against fixed core order. No phylogenetic inference rerun.','limitations':['Core lacks comparable bootstrap support.','Three-gene-derived groups are fixed hypotheses: their recovery in that same tree is circular, not independent validation.','Colors convey inherited three-gene assignments, not independent marker identification.','Literal quartet results count short-edge binary resolutions; RF edge-collapse sensitivity is reported separately.','Bootstrap frequencies conditional on the matching restricted cohort can differ from full reference-rich clade supports.','Unrooted splits are root-independent group separations, not proof of rooted ancestry or nomenclatural validation.','The seven marker combinations overlap in loci and are not seven independent votes.']})
assert all(sha(Path(p))==h for p,h in INPUTS.items())
shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
print('DONE',json.dumps(summary),flush=True)
