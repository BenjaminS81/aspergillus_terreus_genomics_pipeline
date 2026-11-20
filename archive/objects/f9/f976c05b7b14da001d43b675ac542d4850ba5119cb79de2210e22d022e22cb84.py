"""Independent Biopython and sequence/figure checks for the new release."""
from pathlib import Path
import sys,json,hashlib,itertools,random,shutil,re
ROOT=Path(__file__).resolve().parent
from circos_review_common import constrain_process,fasta
limits=constrain_process()
sys.path.append(str(ROOT/'comparative_visualization/2026-11-13/vendor'))
from Bio import Phylo
import fitz
OUT=ROOT/'stakeholder_release/2026-11-20_v4_clade_congruence_cgmlst'
def read(p):return json.loads(p.read_text())
man=read(OUT/'data/analysis_manifest.json');ann=read(OUT/'data/core_tip_assignments.json');names=sorted(a['tip'] for a in ann);alltips=set(names)
def splits(tree):
 result=set()
 for c in tree.find_clades():
  a=frozenset(t.name for t in c.get_terminals());b=frozenset(alltips-a)
  if 1<len(a)<65:result.add(min((a,b),key=lambda s:(len(s),tuple(sorted(s)))))
 return result
trees={k:Phylo.read(OUT/'data'/(k+'.matched66.treefile'),'newick') for k in man['tree_keys']};sets={k:splits(t) for k,t in trees.items()}
summary=read(OUT/'data/global_comparisons.json')
for row in summary:
 a=sets['Core'];b=sets[row['tree']];assert abs(len(a^b)/(len(a)+len(b))-row['normalized_RF'])<1e-12
rng=random.Random(20261120);quartets=[sorted(rng.sample(names,4)) for _ in range(100)]
for key,t in trees.items():
 for q in quartets:
  possible=set()
  for a in sets[key]:
   z=set(q)&a
   if len(z)==2:
    pair=z if q[0] in z else set(q)-z;possible.add(next(i for i in (1,2,3) if q[i] in pair)-1)
  assert len(possible)<=1
  # Independently derived four-point topology from Biopython unit-edge distances.
  old=[c.branch_length for c in t.find_clades()]
  for c in t.find_clades():c.branch_length=1
  a,b,c,d=q;s=[t.distance(a,b)+t.distance(c,d),t.distance(a,c)+t.distance(b,d),t.distance(a,d)+t.distance(b,c)]
  code=s.index(min(s)) if s.count(min(s))==1 else None
  for node,v in zip(t.find_clades(),old):node.branch_length=v
  assert possible==({code} if code is not None else set()),(key,q,s,possible)
original=Phylo.read(ROOT/'stakeholder_release/2026-11-17_v1/data/core_genome_original.treefile','newick')
for t in original.get_terminals():t.name='LOCAL_'+t.name.removesuffix('_transcripts')
assert splits(original)==sets['Core']
delta=max(abs(original.distance(a,b)-trees['Core'].distance(a,b)) for a,b in itertools.combinations(names,2));assert delta<1e-9
p=ROOT/'comparative_visualization/2026-11-11/core66_combinations';current=dict(fasta(ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene/inputs/complete230.fasta'));ranges={'BenA':(0,515),'CaM':(515,1062),'RPB2':(1062,2076)}
for k in man['tree_keys'][1:-1]:
 seq=dict(fasta(p/(k+'.fasta')));assert set(seq)==alltips
 for t in names:assert seq[t]==''.join(current[t][slice(*ranges[g])] for g in k.split('_'))
figures=[]
for pdf in sorted((OUT/'figures').glob('*.pdf')):
 with fitz.open(pdf) as doc:
  assert len(doc)==1;txt=doc[0].get_text();d={'file':pdf.name,'pages':1}
  if 'tanglegram' in pdf.name:
   # Preserve both full labels, including original source names and actual type cultures.
   counts={a['tip']:sum(line.strip()==a['display_label'] for line in txt.splitlines()) for a in ann};assert all(v==2 for v in counts.values()),counts
   assert txt.count('TYPE:')==6;d['full_matching_labels']=132;d['type_labels']=6
  figures.append(d)
for p,h in man['input_hashes'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
# Verify gallery targets and supporting files exist.
doc=(OUT/'review.html').read_text(encoding='utf-8')
for target in re.findall(r'(?:href|src)="([^"]+)"',doc):
 if not target.startswith('https:'):assert (OUT/target).is_file(),target
result={'state':'passed','independent_RF_check':'Biopython-derived splits match every reported RF','quartet_cross_check':'100 random quartets per tree: split-based resolution agrees with independent Biopython four-point distances','original_core_topology_preserved':True,'original_core_max_pairwise_distance_change':delta,'all_seven_marker_alignment_inputs_identical_to_current_matching_sequence_blocks':True,'all_hashed_source_files_unchanged':True,'figure_checks':figures,'gallery_links_valid':True,'limits':limits}
(OUT/'qa/independent_checks.json').write_text(json.dumps(result,indent=2));shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
print(json.dumps({k:v for k,v in result.items() if k!='figure_checks'},indent=2))
