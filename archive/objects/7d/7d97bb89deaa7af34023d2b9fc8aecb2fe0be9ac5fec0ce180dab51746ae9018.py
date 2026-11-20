"""Distinguish exact/subset contig sets from genuine sequence discrepancies."""
from pathlib import Path
import sys,json,hashlib,collections
OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
constrain_process()
def load(name):return json.loads((OUT/'data'/name).read_text())
def save(name,v):(OUT/'data'/name).write_text(json.dumps(v,indent=2))
manifest=load('analysis_eligibility.json');ann=load('annotation_candidates.json');cache={}
def contigs(path):
 path=str(path)
 if path in cache:return cache[path]
 rows=[];name=None;hh=hashlib.sha256();length=0
 with Path(path).open() as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:rows.append({'id':name,'sha256':hh.hexdigest(),'length':length})
    name=line[1:].split()[0];hh=hashlib.sha256();length=0
   else:
    s=line.strip().upper();hh.update(s.encode());length+=len(s)
 if name is not None:rows.append({'id':name,'sha256':hh.hexdigest(),'length':length})
 cache[path]=rows;return rows
canonical={m['sample']:contigs(m['assembly']) for m in manifest}
relationships=[]
for i,a in enumerate(ann):
 if not a['scaffolds_file']:continue
 ar=contigs(a['scaffolds_file']);ac=collections.Counter(r['sha256'] for r in ar)
 for s in a['samples_from_path']:
  cr=canonical[s];cc=collections.Counter(r['sha256'] for r in cr)
  relation='exact_contig_multiset' if ac==cc else 'unchanged_contig_subset' if ac<=cc else 'sequence_or_contig_boundaries_differ'
  byhash=collections.defaultdict(list)
  for r in cr:byhash[r['sha256']].append(r['id'])
  relationships.append({'sample':s,'protein_file':a['protein_file'],'annotation_assembly':a['scaffolds_file'],'relation':relation,'canonical_contigs':len(cr),'annotation_contigs':len(ar),'canonical_bp':sum(r['length'] for r in cr),'annotation_bp':sum(r['length'] for r in ar),'annotation_contigs_not_matched':sum((ac-cc).values()),'canonical_contigs_not_retained':sum((cc-ac).values()),'contig_mapping':[{'annotation_contig':r['id'],'canonical_candidates':byhash[r['sha256']]} for r in ar],'status':'input relationship only; CDS model and excluded-contig effects require review'})
 if (i+1)%25==0:print('Contig relationships',i+1,flush=True)
save('annotation_assembly_relationships.json',relationships)
counts=collections.Counter(r['relation'] for r in relationships)
save('annotation_assembly_relationship_summary.json',{'relationships':dict(counts),'scope':'case-normalised per-contig sequence comparison; orientation/boundary differences are not silently normalized','biological_gene_model_validation':'not performed'})
print(json.dumps(dict(counts)),flush=True)
