"""Build a source-only publication staging tree; never mutate scientific inputs/jobs."""
from pathlib import Path
import os,sys,json,re,csv,hashlib,ast,collections,datetime,shutil
ROOT=Path(__file__).resolve().parent
from circos_review_common import constrain_process
LIMITS=constrain_process()
PROJECT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
OUT=ROOT/'repository_publication/aspergillus_terreus_genomics_pipeline'
OUT.mkdir(parents=True,exist_ok=True)
DATE='2026-11-20'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def digest(b):return hashlib.sha256(b).hexdigest()
def write(rel,text):
 p=OUT/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text,encoding='utf-8',newline='\n')
def js(rel,value):write(rel,json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def tsv(rel,rows,fields):
 p=OUT/rel;p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
  for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})
def location(v):
 s=str(v).replace('\\\\','/').replace('\\','/')
 replacements=[(str(PROJECT).replace('\\','/'),'PROJECT'),('/mnt/d/Illumina/David_Eisele/terreus_genotyping','PROJECT'),(str(ROOT).replace('\\','/'),'WORKSPACE'),('/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics','WORKSPACE')]
 for old,new in replacements:s=s.replace(old,new)
 return s
def clean(v):
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items() if k not in {'sequence','seq','first_headers','raw_xml','xml','raw','fasta','alignment'} and not re.search(r'(?:_sequence|_Seq_trimmed)$',k)}
 if isinstance(v,list):return [clean(x) for x in v]
 return location(v) if isinstance(v,str) else v
SKIP={'vendor','node_modules','busco_downloads','__pycache__','.git','.snakemake','.venv','envs','repository_publication','tmp'}
EXT={'.py','.sh','.bash','.R','.r','.ps1','.nf','.smk'}
versions=[];exclusions=[];syntax=[];objects={}
def category(name):
 n=name.lower()
 if any(x in n for x in ('lov','circos','plain_terreus')):return 'biosynthetic_clusters'
 if any(x in n for x in ('busco','qc','assembly_metrics','transcript_metrics')):return 'quality_control'
 if any(x in n for x in ('ncbi','ena','public','fetch','download','citation','bibliograph','sources')):return 'provenance_retrieval'
 if any(x in n for x in ('plot','render','figure','release','visual','legend')):return 'figures_reporting'
 if any(x in n for x in ('tree','phylo','monophy','clade','core66','congru','combined','concaten','type_')):return 'phylogenetics'
 if any(x in n for x in ('ortholog','orthogroup','marker','locus','reference','taxonom','genomic')):return 'taxonomy_orthology'
 if any(x in n for x in ('gate','checkpoint','capped','resource','pause','resume')):return 'execution_gates'
 return 'inventory_utilities'
def archive(p,root,origin,content=None,transformation='none'):
 raw=p.read_bytes();b=raw if content is None else content.encode('utf-8');sha=digest(b);rel=p.relative_to(root).as_posix()
 # Content-addressed storage preserves every distinct version without repeated copies.
 stored=objects.get(sha)
 if stored is None:
  suffix=p.suffix if p.suffix in EXT|{'.yaml','.yml'} and transformation=='none' else '.txt'
  stored=f'archive/objects/{sha[:2]}/{sha}{suffix}';dest=OUT/stored;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b);objects[sha]=stored
 row={'origin':origin,'original_relative_path':rel,'stored_path':stored,'original_sha256':digest(raw),'stored_sha256':sha,'bytes':len(b),'source_mtime_utc':datetime.datetime.fromtimestamp(p.stat().st_mtime,datetime.timezone.utc).isoformat(),'snapshot_version':next((x for x in p.parts if re.match(r'^2026-\d\d-\d\d',x)),'undated_original'),'category':category(p.name),'transformation':transformation,'latest_workspace_entry':origin=='workspace' and p.parent==ROOT}
 versions.append(row)
 if p.suffix=='.py' and transformation=='none':
  try:ast.parse(b.decode('utf-8-sig'));state='parses'
  except (SyntaxError,UnicodeError) as e:state=f'{type(e).__name__}: {str(e).splitlines()[0]}'
  syntax.append({'source':origin+'/'+rel,'stored_path':stored,'status':state})
 if row['latest_workspace_entry']:
  target=OUT/'scripts/current'/p.name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(b)
 return row
for folder,dirs,files in os.walk(ROOT):
 dirs[:]=[d for d in dirs if d not in SKIP]
 for name in files:
  p=Path(folder)/name
  if p.suffix in EXT and p.name not in {'prepare_git_publication.py'}:archive(p,ROOT,'workspace')
# Original hand-maintained commands, environments and notebook code. Never publish notebook outputs.
for sub in ('workflow','scr','terreus_qc_check/scripts','config/env','docs/steps'):
 for p in sorted((PROJECT/sub).rglob('*')):
  if not p.is_file():continue
  if p.suffix in {'.odt','.odg','.pdf','.png','.json','.tsv','.log'}:exclusions.append({'origin':'project','path':p.relative_to(PROJECT).as_posix(),'reason':'non-source document, result or metadata; inventory only','sha256':digest(p.read_bytes())});continue
  if p.suffix in EXT|{'.yaml','.yml','.txt'} or not p.suffix:archive(p,PROJECT,'project')
for p in sorted((PROJECT/'notebooks').glob('*.ipynb')):
 nb=read(p);chunks=[]
 for i,c in enumerate(nb.get('cells',[])):
  if c.get('cell_type')=='code':chunks.append(f'# Original notebook code cell {i}\n'+''.join(c.get('source',[])))
 archive(p,PROJECT,'project',content='\n\n'.join(chunks)+'\n',transformation='code cells only; outputs, attachments, markdown and execution metadata omitted')
# A mixed lab notebook contains results and narrative: preserve its source hash and export fenced command blocks only.
nb=PROJECT/'docs/lab-notebook.md'
if nb.exists():
 text=nb.read_text(encoding='utf-8');blocks=re.findall(r'^```([^\n]*)\n(.*?)^```\s*$',text,re.M|re.S)
 chosen=[f'# Command block {i}; original fence: {lang}\n{code}' for i,(lang,code) in enumerate(blocks) if lang.strip().lower() in {'bash','sh','shell','python','r','powershell'}]
 if chosen:archive(nb,PROJECT,'project',content='\n\n'.join(chosen),transformation='language-tagged executable code fences only; mixed narrative/results omitted')
 exclusions.append({'origin':'project','path':'docs/lab-notebook.md','reason':'Full mixed notebook omitted because it contains results; language-tagged code fences exported','sha256':digest(nb.read_bytes())})
cfg=PROJECT/'config/config.yaml';archive(cfg,PROJECT,'project')
js('manifests/script_versions.json',versions);js('manifests/archive_exclusions.json',exclusions);js('manifests/historical_python_syntax.json',syntax)
print('SOURCE',len(versions),'records',len(objects),'distinct objects',flush=True)
# Scientific inventories: explicitly selected fields only. No raw sequence strings.
T=ROOT/'taxonomy_reference_build/2026-11-10';LIB=T/'reference_library';A=T/'analysis';AUDIT=ROOT/'review_reports/2026-11-20_public_genome_audit_v1'
refs=read(LIB/'validated_reference_rows.json');gb=read(LIB/'genbank_records.json');analysisrefs={x['reference_id']:x for x in read(A/'tree_reference_metadata.json')}
from terreus_tree_helpers import parse_newick
treepaths={'core66':ROOT/'stakeholder_release/2026-11-17_v1/data/core_genome_original.treefile','complete230':ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene/runs/complete230.treefile','full283':A/'trees/combined_all_strains.treefile','type108':A/'trees/type_focused_sensitivity.treefile'}
tipsets={k:{n.name for n in parse_newick(p.read_text()).leaves()} for k,p in treepaths.items()}
refrows=[]
for r in refs:
 ar=analysisrefs.get(r['reference_id'],{});refrows.append({k:r.get(k,'') for k in ['reference_id','species','strain','type_flag','section','series','BenAGB','CaMGB','RPB2GB']}|{'analysis_species':ar.get('analysis_species',''),'in_full283':r['reference_id'] in tipsets['full283'],'in_complete230':r['reference_id'] in tipsets['complete230'],'in_type108':r['reference_id'] in tipsets['type108']})
tsv('inventory/reference_strains.tsv',refrows,list(refrows[0]))
gbrows=[{'accession':k,'accession_version':r.get('accession-version',''),'organism':r.get('organism',''),'definition':r.get('definition',''),'length':r.get('length',''),'status':'sequence fetched; payload excluded from repository','url':'https://www.ncbi.nlm.nih.gov/nuccore/'+r.get('accession-version',k)} for k,r in sorted(gb.items())]
tsv('inventory/fetched_marker_accessions.tsv',gbrows,list(gbrows[0]))
calls=read(ROOT/'review_reports/2026-11-20_public_genome_audit_v2/actual_tree_inclusion_all_sources.json');canon=read(ROOT/'gate_preparation/2026-11-18_v1/data/canonical_assemblies.json');ch={s:r for r in canon for s in r['samples']}
assign={x['sample']:x for x in read(ROOT/'stakeholder_release/2026-11-20_v3_original_taxon_labels/data/core_tip_assignments.json')}
ann=read(ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene/data/new_tree_annotations.json')
genomes=[]
for r in calls:
 s=r['sample'];a=assign.get(s,{});c=ch.get(s,{});genomes.append({'sample':s,'source_taxon':a.get('reported_species',ann.get('LOCAL_'+s,{}).get('taxon','unknown')),'type_derived_assignment':a.get('assigned_species','not in core annotation'),'assignment_evidence':a.get('confidence','unknown'),'public_data':r.get('public_data'),'provenance':r.get('provenance'),'run_accessions':r.get('run_accessions',[]),'assembly_accessions':r.get('assembly_accessions',[]),'assembly_location':location(r.get('assembly_path','')),'assembly_present_audited':r.get('assembly_present'),'sequence_sha256':c.get('sequence_sha256',''),'biological_culture_equivalence':'not established by sequence identity','core66':r.get('core66'),'complete230':r.get('complete230'),'full283':r.get('full283')})
tsv('inventory/genome_cohort.tsv',genomes,list(genomes[0]))
assemblies=read(AUDIT/'ncbi_assemblies.json');am={x['accession']:x for x in assemblies};cross=read(ROOT/'review_reports/2026-11-20_public_genome_audit_v2/assembly_relevance_crosscheck.json')
js('inventory/public_assembly_crosscheck.json',clean(cross))
rawfiles=read(AUDIT/'local_raw_read_files.json');rawrows=[{'location':location(r['path']),'run_accession':r.get('run_accession',''),'bytes':r.get('bytes'),'observed_mtime_ns':r.get('mtime_ns'),'scope':'file presence inventory; not proof of read integrity or assembly linkage'} for r in rawfiles]
tsv('inventory/local_raw_read_files.tsv',rawrows,list(rawrows[0]))
ena=read(AUDIT/'ena_read_run.json');runfields=['run_accession','study_accession','sample_accession','scientific_name','strain','isolate','library_strategy','library_source','library_layout','instrument_platform','read_count','base_count','fastq_ftp','fastq_md5','fastq_bytes'];localruns={r.get('run_accession') for r in rawfiles}
enarows=[{k:r.get(k,'') for k in runfields}|{'local_file_observed':r.get('run_accession') in localruns} for r in ena]
tsv('inventory/ena_runs_queried.tsv',enarows,runfields+['local_file_observed'])
added=read(ROOT/'review_reports/2026-11-20_public_genome_audit_v3_fetch/download_manifest.json')
js('inventory/public_additions.json',clean(added))
downloadrows=[]
for r in read(T/'download_manifest.json'):downloadrows.append(clean(r)|{'manifest':'taxonomy_reference_build/2026-11-10/download_manifest.json'})
js('inventory/publication_downloads.json',downloadrows)
js('inventory/publications.json',read(LIB/'marker_publications_enriched.json'))
# Assembly and transcript source asset inventories retain paths/checksums/status, never payloads.
for src,dst in [(A/'assembly_asset_inventory.json','assembly_assets.json'),(A/'assembly_variant_inventory.json','assembly_variants.json'),(ROOT/'gate_preparation/2026-11-18_v1/data/analysis_eligibility.json','analysis_eligibility.json')]:
 if src.exists():js('inventory/'+dst,clean(read(src)))
# Taxon evidence catalog: preserve source names, synonyms and provisional names; do not silently synonymize.
taxa=collections.defaultdict(lambda:collections.defaultdict(set))
def taxon(name,role,evidence):
 if not name:return
 name=str(name).replace('_',' ').strip();name=re.sub(r'^A\.\s*','Aspergillus ',name);name=re.sub(r'\s+',' ',name)
 m=re.match(r'^([A-Z][a-z]+) ([a-z][a-z-]+)(?: (var\.|subsp\.) ([a-z-]+))?',name)
 if not m:return
 normalized=' '.join(x for x in m.groups() if x);taxa[normalized][role].add(str(evidence))
for r in refs:
 taxon(r.get('species'),'reference_library',r['reference_id'])
 for label in ('full283','complete230','type108'):
  if r['reference_id'] in tipsets[label]:taxon(r.get('species'),'used_'+label,r['reference_id'])
for r in analysisrefs.values():
 for label in ('full283','complete230','type108'):
  if r['reference_id'] in tipsets[label]:taxon(r.get('analysis_species'),'analysis_name_'+label,r['reference_id'])
for k,r in gb.items():taxon(r.get('organism'),'marker_sequence_fetched',k)
for r in read(LIB/'accepted_species.json'):taxon(r.get('ID'),'nomenclature_catalog_only',r.get('MYCOBANK',''))
for r in genomes:
 taxon(r['source_taxon'],'genome_source_name',r['sample'])
 if r['core66']:
  taxon(r['source_taxon'],'used_core_source_name',r['sample']);taxon(r['type_derived_assignment'],'provisional_core_assignment',r['sample'])
for r in assemblies:taxon(r.get('organism',{}).get('organism_name'),'assembly_metadata_queried',r['accession'])
for r in ena:taxon(r.get('scientific_name'),'read_metadata_queried',r['run_accession'])
for r in added['assemblies']:
 organism=am.get(r['accession'],{}).get('organism',{}).get('organism_name')
 if organism:taxon(organism,'genome_download_verified',r['accession'])
taxrows=[]
for name,evidence in sorted(taxa.items()):
 taxrows.append({'taxon_name':name,'genus':name.split()[0],'sequence_fetched':bool(evidence.get('marker_sequence_fetched') or evidence.get('genome_download_verified')),'used_core':bool(evidence.get('used_core_source_name') or evidence.get('provisional_core_assignment')),'used_complete230':bool(evidence.get('used_complete230') or evidence.get('analysis_name_complete230')),'used_full283':bool(evidence.get('used_full283') or evidence.get('analysis_name_full283')),'roles':'; '.join(sorted(evidence)),'evidence':{k:sorted(v) for k,v in sorted(evidence.items())}})
js('inventory/taxon_evidence.json',taxrows)
tsv('inventory/taxa.tsv',[{k:v for k,v in r.items() if k!='evidence'} for r in taxrows],['taxon_name','genus','sequence_fetched','used_core','used_complete230','used_full283','roles'])
genera=sorted({r['genus'] for r in taxrows});write('inventory/genera.txt','\n'.join(genera)+'\n')
write('inventory/species_and_infraspecific_names.txt','\n'.join(r['taxon_name'] for r in taxrows)+'\n')
summary={'inventory_date':DATE,'script_source_records':len(versions),'distinct_archived_source_objects':len(objects),'current_workspace_scripts':sum(r['latest_workspace_entry'] for r in versions),'python_archive_syntax_exceptions':sum(r['status']!='parses' for r in syntax),'reference_library_strains':len(refs),'fetched_marker_accessions':len(gb),'genome_cohort_labels':len(genomes),'canonical_assembly_sequences':len(canon),'local_raw_read_files':len(rawrows),'ena_runs_queried':len(ena),'taxon_names_all_scopes':len(taxrows),'genera':genera,'taxon_names_with_fetched_sequences':sum(r['sequence_fetched'] for r in taxrows),'taxon_names_used_complete230_including_source_and_analysis_names':sum(r['used_complete230'] for r in taxrows),'new_assemblies_verified':len(added['assemblies']),'limits':LIMITS,'scope':'Metadata and code only. Catalog-only and queried-only taxa are not represented as downloaded sequences. Available snapshots are preserved; unrecorded historical edits cannot be reconstructed.'}
js('manifests/publication_summary.json',summary)
js('manifests/tree_versions.json',[{'version':k,'tips':len(tipsets[k]),'source':location(p),'sha256':digest(p.read_bytes()),'payload_in_repository':False} for k,p in treepaths.items()])
snapshots=[]
for p in sorted(ROOT.glob('*/*')):
 if p.is_dir() and re.match(r'^2026-\d\d-\d\d',p.name):snapshots.append({'path':p.relative_to(ROOT).as_posix(),'date_label':p.name,'archived_script_records':sum(r['origin']=='workspace' and r['original_relative_path'].startswith(p.relative_to(ROOT).as_posix()+'/') for r in versions),'payloads_excluded':True})
js('manifests/analysis_versions.json',snapshots)
# Snapshot gate status is dated evidence, never a claim that an active job is complete.
for src,name in [(ROOT/'gate_preparation/2026-11-18_v1/data/outstanding_gate_actions.json','initial_gate_snapshot.json'),(ROOT/'stakeholder_release/2026-11-20_v1_type_monophyly/data/evodevo_gates.json','evodevo_gates.json'),(ROOT/'stakeholder_release/2026-11-20_v4_clade_congruence_cgmlst/data/cgmlst_gates.json','cgmlst_gates.json'),(ROOT/'gate_execution/2026-11-18_v1/status.json','execution_status_snapshot.json')]:
 if src.exists():js('manifests/'+name,clean(read(src)))
write('docs/TAXA.md','# Taxon coverage\n\nThis is an evidence catalog, not a taxonomic revision. Source and analysis names are retained separately, including synonyms and provisional assignments. Catalog-only names are not claimed to have been sequenced or downloaded.\n\nGenera: '+', '.join(genera)+'.\n\nSee [all names](../inventory/species_and_infraspecific_names.txt), [scope table](../inventory/taxa.tsv), and [record-level evidence](../inventory/taxon_evidence.json).\n\n| Name | Sequence fetched | Used in current complete three-gene tree |\n|---|---|---|\n'+''.join(f"| {r['taxon_name']} | {r['sequence_fetched']} | {r['used_complete230']} |\n" for r in taxrows))
write('docs/SCRIPT_CATALOG.md','# Script catalog\n\nCurrent means the observed top-level workspace version, not a certification of portability or biological validity. All distinct historical snapshots are indexed by content hash. Run archival workflows only in a separately restored analysis workspace after reviewing their paths and resource limits.\n\n| Category | Current script | Source hash |\n|---|---|---|\n'+''.join(f"| {r['category']} | [{Path(r['original_relative_path']).name}](../scripts/current/{Path(r['original_relative_path']).name}) | `{r['stored_sha256'][:12]}` |\n" for r in sorted(versions,key=lambda x:(x['category'],x['original_relative_path'])) if r['latest_workspace_entry']))
print(json.dumps(summary,indent=2),flush=True)
