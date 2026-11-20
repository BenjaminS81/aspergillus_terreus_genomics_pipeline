"""Resolve boundary-only assembly differences using unique exact non-N contig bodies.
Never infer a relationship from a similar sequence or a shared strain name alone.
"""
from pathlib import Path
import hashlib,collections,json
from worker import OUT,PREP,lp,load,save,sha
from validate_annotations import seqdict
def digest(s):return hashlib.sha256(s.encode()).hexdigest()
def main():
 proposals=load(PREP/'data/annotation_release_proposals.json');candidates=load(PREP/'data/annotation_candidates.json');manifest={m['sample']:m for m in load(PREP/'data/analysis_eligibility.json')};evidence=[]
 for p in proposals:
  if p['preferred_candidate']:continue
  s=p['sample'];canonical_path=lp(manifest[s]['assembly']);canonical=seqdict(canonical_path)
  assert sha(canonical_path)==manifest[s]['assembly_file_sha256']
  idx=collections.defaultdict(list)
  for cid,seq in canonical.items():
   body=seq.strip('N')
   if body:idx[digest(body)].append((cid,len(seq)-len(seq.lstrip('N'))))
  options=[a for a in candidates if s in a['samples_from_path'] and a['scaffolds_file'] and a['gff'] and not a['unmapped_proteins']]
  options.sort(key=lambda a:(a['stage']=='annotate',a['mtime_ns']),reverse=True)
  for a in options:
   genome=seqdict(lp(a['scaffolds_file']));mapping={};bad=[]
   for cid,seq in genome.items():
    body=seq.strip('N');hits=idx.get(digest(body),[]) if body else []
    if len(hits)!=1:bad.append({'contig':cid,'matches':len(hits)});continue
    target,left=hits[0];offset=left-(len(seq)-len(seq.lstrip('N')))
    mapping[cid]={'canonical_contig':target,'offset':offset}
   # Every annotated interval must map exactly, including strand-independent sequence.
   gff_bad=[]
   if not bad:
    with lp(a['gff']).open() as f:
     for n,line in enumerate(f,1):
      if line.startswith('#'):continue
      ff=line.rstrip().split('\t')
      if len(ff)!=9:continue
      cid=ff[0]
      if cid not in mapping:gff_bad.append(n);continue
      start,end=int(ff[3])-1,int(ff[4]);m=mapping[cid];cs,ce=start+m['offset'],end+m['offset']
      if cs<0 or ce>len(canonical[m['canonical_contig']]) or genome[cid][start:end]!=canonical[m['canonical_contig']][cs:ce]:gff_bad.append(n)
   report={'sample':s,'candidate':a['protein_file'],'annotation_genome_sha256':sha(lp(a['scaffolds_file'])),'gff_sha256':sha(lp(a['gff'])),'canonical_sha256':sha(canonical_path),'unmatched_or_ambiguous_contigs':bad,'invalid_feature_rows':gff_bad,'mapping':mapping,'status':'EXACT_BOUNDARY_MAPPING_PASS' if not bad and not gff_bad else 'UNRESOLVED'}
   evidence.append(report)
   if report['status']=='EXACT_BOUNDARY_MAPPING_PASS':
    p.update({'preferred_candidate':a['protein_file'],'assembly_relationship':'unique_exact_contig_bodies_with_terminal_N_offsets','status':'boundary mapping passed; now requires CDS translation validation'});break
  print(s,'resolved' if p['preferred_candidate'] else 'unresolved',flush=True)
 save(OUT/'data/resolved_annotation_proposals.json',proposals);save(OUT/'reports/annotation_boundary_resolution.json',evidence)
if __name__=='__main__':main()
