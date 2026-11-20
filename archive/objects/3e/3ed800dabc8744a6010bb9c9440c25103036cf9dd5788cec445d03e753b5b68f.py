from pathlib import Path
import json,hashlib,shutil,csv,collections
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene'
OLD=ROOT/'taxonomy_reference_build/2026-11-10/analysis'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fasta(p):
 d={}
 for line in p.read_text().splitlines():
  if line.startswith('>'):k=line[1:].split()[0];assert k not in d;d[k]=''
  elif line.strip():d[k]+=line.strip()
 return d
for name in ['inputs','data','runs','figures','scripts','qa']:(OUT/name).mkdir(parents=True,exist_ok=True)
seq=fasta(OLD/'tree_inputs/combined_all_strains.fasta');occ=read(OLD/'tree_inputs/locus_occupancy.json')
refs={r['reference_id']:r for r in read(OLD/'tree_reference_metadata.json')}
keep={r['id'] for r in occ if r['available_loci']==3}
assert len(keep)==230 and keep<=set(seq)
genes=[('BenA',0,515),('CaM',515,1062),('RPB2',1062,2076)]
rows=[]
for r in occ:
 k=r['id'];row={'tip':k,'retained':k in keep,'source':'genome' if k.startswith('LOCAL_') else 'reference','available_loci':r['available_loci']}
 for g,a,b in genes:
  n=sum(c in 'ACGT' for c in seq[k][a:b].upper());assert n==r[g+'_ACGT'];row[g+'_ACGT']=n
  if k in keep:assert n>=120 and set(seq[k][a:b])!={'?'}
 row['reference_type']=bool(refs.get(k,{}).get('type_flag'));row['reference_species']=refs.get(k,{}).get('analysis_species','')
 rows.append(row)
(OUT/'inputs/complete230.fasta').write_text(''.join('>'+k+'\n'+seq[k]+'\n' for k in sorted(keep)))
shutil.copy2(OLD/'tree_inputs/partitions.nex',OUT/'inputs/partitions.nex')
with (OUT/'data/selection_audit.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
full=OUT/'runs/full283_existing';full.mkdir(exist_ok=True)
original=OLD/'trees/combined_all_strains'
plot_tree=ROOT/'comparative_visualization/2026-11-11/sources/input_trees/combined_all_strains.treefile'
assert sha(Path(str(original)+'.treefile'))==sha(plot_tree)
from terreus_tree_helpers import parse_newick
assert {t.name for t in parse_newick(plot_tree.read_text()).leaves()}==set(seq)
boot=Path(str(original)+'.ufboot');n=0
for line in boot.open():
 if not line.strip():continue
 assert {t.name for t in parse_newick(line).leaves()}==set(seq);n+=1
assert n==1000
report=Path(str(original)+'.iqtree').read_text();assert '1000 replicates' in report and 'SH-aLRT' in report
hashes={}
for ext in ['.treefile','.contree','.ufboot','.iqtree','.log','.best_scheme.nex','.splits.nex']:
 src=Path(str(original)+ext);dest=full/src.name;shutil.copy2(src,dest);hashes[str(src)]=sha(src)
for name in ['01_combined_three_gene_vertical.pdf','02_combined_three_gene_circular.pdf','04_consolidated_colour_legend.pdf']:
 shutil.copy2(ROOT/'stakeholder_release/2026-11-17_v3_taxonomy_review'/name,OUT/'figures'/('full283_'+name))
design={'full_tree':'existing completed IQ-TREE 2.4.0 analysis; verified, not rerun','full_tips':283,'verified_full_bootstraps':n,'complete_tips':len(keep),'genomes':sum(k.startswith('LOCAL_') for k in keep),'references':sum(not k.startswith('LOCAL_') for k in keep),'complete_reference_types':sum(bool(r['type_flag']) for k,r in refs.items() if k in keep),'excluded_reference_types':[{ 'tip':k,'species':r['analysis_species']} for k,r in refs.items() if r['type_flag'] and k not in keep],'unrepresented_reference_names':sorted({r['analysis_species'] for r in refs.values()}-{r['analysis_species'] for k,r in refs.items() if k in keep}),'selection':'All existing entries with three represented marker regions; no taxon downsampling or imputation. Existing curated locus assignments retained; complete occupancy is not new proof of single-copy orthology.','sites':2076,'outgroup':'REF_436_Aspergillus_flavipes_NRRL302','source_hashes':hashes,'alignment_sha256':sha(OUT/'inputs/complete230.fasta'),'partitions_sha256':sha(OUT/'inputs/partitions.nex')}
(OUT/'data/design.json').write_text(json.dumps(design,indent=2));print(json.dumps({k:v for k,v in design.items() if k!='source_hashes'},indent=2))
