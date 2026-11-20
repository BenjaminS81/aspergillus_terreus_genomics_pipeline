"""Repair legacy counting units, preserving explicit input/database uncertainty."""
from pathlib import Path
import json,csv,re,collections,math
from worker import OUT,PREP,lp,load,save,sha
PROJECT=Path('/mnt/d/Illumina/David_Eisele/terreus_genotyping')
FAMILY=re.compile(r'(?<![A-Za-z])(?:AA|CBM|CE|GH|GT|PL)\d+')
DOMAIN=re.compile(r'((?:AA|CBM|CE|GH|GT|PL)\d+(?:_\d+)?)\((\d+)-(\d+)\)')
def families(s):return set(FAMILY.findall(s or ''))
def tests():
 assert families('CBM87(97-314)+CE18(315-682)')=={'CBM87','CE18'}
 assert families('AA3_2+AA3_2')=={'AA3'}
 assert len(DOMAIN.findall('CBM87(97-314)+CE18(315-682)'))==2
 save(OUT/'reports/functional_parser_tests.json',{'multidomain_split':'PASS','gene_family_deduplication':'PASS','coordinate_domain_units':'PASS'})
def readtable(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
 tests();ann=load(OUT/'reports/annotation_release_summary.json');rows=[];events=[];eff=[]
 for a in ann:
  s=a['sample'];validfile=OUT/'runs/annotation_release'/s/'validation.json'
  if not validfile.exists():rows.append({'sample':s,'status':'annotation_unresolved','counts':None});continue
  v=load(validfile);primary={r['protein_id'] for r in v['model_rows'] if r['primary_selected']};allids={r['protein_id'] for r in v['model_rows']}
  p=PROJECT/'results/0_8_dbcan'/s/'overview.tsv'
  if not p.exists():rows.append({'sample':s,'status':'no_sample_matched_legacy_dbCAN','counts':None});continue
  rr=readtable(p);foreign={r['Gene ID'] for r in rr}-allids
  if foreign:rows.append({'sample':s,'status':'foreign_protein_ids','foreign_count':len(foreign),'counts':None});continue
  counters=collections.Counter();gene_family=set();domain_counter=collections.Counter();gene_supported=set();domain_occurrences=set()
  for r in rr:
   pid=r['Gene ID']
   if pid not in primary:continue
   calls={tool:families(r.get(tool,'')) for tool in ['dbCAN_hmm','dbCAN_sub','DIAMOND']}
   if int(r.get('#ofTools','0'))>=2:gene_supported.add(pid)
   for family in set.union(*calls.values()):
    methods=[tool for tool,ff in calls.items() if family in ff]
    if len(methods)>=2:gene_family.add((pid,family))
    events.append({'sample':s,'protein':pid,'family':family,'methods':methods,'family_supported_by_at_least_two_methods':len(methods)>=2})
   for family,start,end in DOMAIN.findall(r.get('dbCAN_hmm','')):
    domain_occurrences.add((pid,family,int(start),int(end)))
  for pid,family,start,end in domain_occurrences:domain_counter[family]+=1
  for pid,family in gene_family:counters[family]+=1
  rows.append({'sample':s,'status':'reparsed_provisional','counts':dict(counters),'HMM_domain_counts':dict(domain_counter),'genes_supported_by_two_tools':len(gene_supported),'distinct_gene_family_pairs':len(gene_family),'input_output_sha256':sha(p),'annotation_status':v['technical_status'],'historical_input_and_database_lock':'not recovered; current ID compatibility is not proof of identical historical sequence input'})
  ep=PROJECT/'results/0_9_effectorp'/s/'effectorp_res'
  if ep.exists():
   er=readtable(ep);ids={r['# Identifier'].split()[0] for r in er};counts=collections.Counter(r['Prediction'] for r in er if r['# Identifier'].split()[0] in primary)
   eff.append({'sample':s,'counts':dict(counts),'legacy_id_set_equals_current_annotation':ids==allids,'output_sha256':sha(ep),'scope':'whole-proteome predictions only; secretome and transmembrane validation not run'})
 save(OUT/'reports/dbcan_reconciled_counts.json',rows);save(OUT/'data/dbcan_gene_family_evidence.json',events);save(OUT/'reports/effectorp_primary_counts.json',eff)
 complete=[r for r in rows if r['counts'] is not None and r.get('annotation_status')=='PASS']
 fams=sorted({f for r in complete for f in r['counts']});matrix=[[r['counts'].get(f,0) for f in fams] for r in complete]
 save(OUT/'data/dbcan_primary_family_matrix.json',{'samples':[r['sample'] for r in complete],'families':fams,'counts':matrix,'unknown_samples':[r['sample'] for r in rows if r not in complete],'status':'provisional legacy evidence; missing rows excluded, never converted to biological zero'})
 # Deterministic average-linkage clustering on binary Jaccard distances. No forced k.
 names=[r['sample'] for r in complete];sets=[{f for f in fams if r['counts'].get(f,0)>0} for r in complete];n=len(names)
 dist={(i,j):(1-len(sets[i]&sets[j])/len(sets[i]|sets[j]) if sets[i]|sets[j] else 0.0) for i in range(n) for j in range(i+1,n)}
 clusters={i:[i] for i in range(n)};links=[];heights={i:0.0 for i in range(n)};newick={i:"'"+names[i].replace("'","''")+"'" for i in range(n)}
 while len(clusters)>1:
  candidates=[];kk=sorted(clusters)
  for ii,i in enumerate(kk):
   for j in kk[ii+1:]:
    dd=sum(dist[tuple(sorted((x,y)))] for x in clusters[i] for y in clusters[j])/(len(clusters[i])*len(clusters[j]));candidates.append((dd,i,j))
  d,i,j=min(candidates);k=n+len(links);members=clusters.pop(i)+clusters.pop(j);clusters[k]=members
  h=d/2
  assert h+1e-10>=max(heights[i],heights[j]);newick[k]=f'({newick[i]}:{max(0,h-heights[i]):.8f},{newick[j]}:{max(0,h-heights[j]):.8f})';heights[k]=h;links.append([i,j,d,len(members)])
 save(OUT/'data/dbcan_jaccard_average_linkage.json',{'samples':names,'linkage':links,'metric':'binary Jaccard','method':'average linkage','Newick_heights':'half the merge distance; branch = parent height minus child height','cluster_count':'not forced','status':'exploratory, legacy prediction provenance unresolved'})
 if n:(OUT/'data/dbcan_jaccard_average_linkage.nwk').write_text(newick[next(iter(clusters))]+';\n')
 save(OUT/'reports/functional_reanalysis_acceptance.json',{'dbcan_status_counts':dict(collections.Counter(r['status'] for r in rows)),'matrix_samples':n,'families':len(fams),'linkage_rows':len(links),'leaf_set_preserved':len({x for cc in clusters.values() for x in cc})==n,'gate_status':'counting/clustering implementation corrected; predictions and biological comparisons remain conditional on provenance/secretome validation'})
 print('Functional counts and clustering regenerated; samples',n,'families',len(fams),'legacy biological limitations retained',flush=True)
if __name__=='__main__':main()
