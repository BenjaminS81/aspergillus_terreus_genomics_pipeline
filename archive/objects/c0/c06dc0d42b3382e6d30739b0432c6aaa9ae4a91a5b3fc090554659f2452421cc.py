"""Complete fetched-accession and auxiliary-reference inventories in staging."""
from pathlib import Path
import json,csv,re,hashlib,shutil
ROOT=Path(__file__).resolve().parent
from circos_review_common import constrain_process
constrain_process()
OUT=ROOT/'repository_publication/aspergillus_terreus_genomics_pipeline';INV=OUT/'inventory';LIB=ROOT/'taxonomy_reference_build/2026-11-10/reference_library'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def tsv(p,rows,fields):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def loc(p):return str(p).replace('\\\\','/').replace('\\','/').replace('D:/Illumina/David_Eisele/terreus_genotyping','PROJECT')
with (INV/'fetched_marker_accessions.tsv').open(encoding='utf-8',newline='') as f:old=list(csv.DictReader(f,delimiter='\t'))
records={r['accession']:r for r in old};historical=read(LIB/'historical_2011/records.json')
for r in records.values():r['retrieval_collections']='curated_reference_library'
for h in historical:
 key=h['primary-accession']
 if key in records:records[key]['retrieval_collections']+='; historical_2011'
 else:records[key]={'accession':key,'accession_version':h.get('accession-version',''),'organism':h.get('organism',''),'definition':h.get('definition',''),'length':len(h.get('sequence','')),'status':'sequence fetched; payload excluded from repository','url':'https://www.ncbi.nlm.nih.gov/nuccore/'+h.get('accession-version',key),'retrieval_collections':'historical_2011'}
rows=[records[k] for k in sorted(records)];tsv(INV/'fetched_marker_accessions.tsv',rows,list(rows[0]))
proteins={}
for name in ('uniprot_lov_genes.json','uniprot_lov_cluster_nine.json'):
 for r in read(ROOT/'comparative_visualization/2026-11-11/sources'/name)['results']:
  proteins[r['primaryAccession']]={'accession':r['primaryAccession'],'organism':r['organism']['scientificName'],'genes':r.get('genes',[]),'length':r['sequence']['length'],'sequence_version':r.get('entryAudit',{}).get('sequenceVersion'),'source_file':name,'url':'https://www.uniprot.org/uniprotkb/'+r['primaryAccession'],'role':'lovastatin cluster protein reference; sequence payload excluded'}
save(INV/'lovastatin_protein_references.json',list(proteins.values()))
taxa={r['taxon_name']:r for r in read(INV/'taxon_evidence.json')}
def add(name,role,acc):
 m=re.match(r'^([A-Z][a-z]+) ([a-z][a-z-]+)(?: (var\.|subsp\.) ([a-z-]+))?',name)
 if not m:return
 name=' '.join(x for x in m.groups() if x)
 row=taxa.setdefault(name,{'taxon_name':name,'genus':name.split()[0],'sequence_fetched':False,'used_core':False,'used_complete230':False,'used_full283':False,'roles':'','evidence':{}})
 row['evidence'].setdefault(role,[])
 if acc not in row['evidence'][role]:row['evidence'][role].append(acc)
 row['sequence_fetched']=True;row['roles']='; '.join(sorted(row['evidence']))
for r in historical:add(r['organism'],'historical_sequence_fetched',r['primary-accession'])
for r in proteins.values():add(r['organism'],'protein_reference_fetched',r['accession'])
taxrows=[taxa[k] for k in sorted(taxa)];save(INV/'taxon_evidence.json',taxrows)
fields=['taxon_name','genus','sequence_fetched','used_core','used_complete230','used_full283','roles'];tsv(INV/'taxa.tsv',[{k:r[k] for k in fields} for r in taxrows],fields)
(INV/'genera.txt').write_text('\n'.join(sorted({r['genus'] for r in taxrows}))+'\n')
(INV/'species_and_infraspecific_names.txt').write_text('\n'.join(taxa)+'\n')
# Ensure deterministic sorting even if an additional historical taxon was inserted.
(INV/'species_and_infraspecific_names.txt').write_text('\n'.join(sorted(taxa))+'\n')
localasm=read(ROOT/'review_reports/2026-11-20_public_genome_audit_v1/local_assembly_files.json')
tsv(INV/'local_public_assembly_files.tsv',[{'location':loc(r['path']),'accession':r.get('accession',''),'bytes':r.get('bytes'),'observed_mtime_ns':r.get('mtime_ns'),'nonempty':r.get('nonempty')} for r in localasm],['location','accession','bytes','observed_mtime_ns','nonempty'])
with (INV/'reference_strains.tsv').open(encoding='utf-8',newline='') as f:refs=list(csv.DictReader(f,delimiter='\t'))
histmeta=read(LIB/'historical_2011/merged_strain_metadata.json')
save(INV/'historical_reference_strains.json',[{k:v for k,v in r.items() if 'sequence' not in k.lower() and 'seq_trimmed' not in k.lower()} for r in histmeta])
extra=read(ROOT/'comparative_visualization/2026-11-11/sources/download_manifest.json')
downloads=read(INV/'publication_downloads.json')
for r in extra:downloads.append(r|{'manifest':'comparative_visualization/2026-11-11/sources/download_manifest.json'})
save(INV/'publication_downloads.json',downloads)
summary=read(OUT/'manifests/publication_summary.json');summary.update(fetched_marker_accessions=len(records),historical_marker_accessions=len(historical),additional_historical_accessions=len(records)-len(old),lovastatin_protein_reference_accessions=len(proteins),local_public_assembly_files=len(localasm),taxon_names_all_scopes=len(taxrows),taxon_names_with_fetched_sequences=sum(r['sequence_fetched'] for r in taxrows))
save(OUT/'manifests/publication_summary.json',summary)
doc=(OUT/'docs/TAXA.md').read_text(encoding='utf-8');prefix=doc.split('| Name |')[0]
(OUT/'docs/TAXA.md').write_text(prefix+'| Name | Sequence fetched | Used in current complete three-gene tree |\n|---|---|---|\n'+''.join(f"| {r['taxon_name']} | {r['sequence_fetched']} | {r['used_complete230']} |\n" for r in taxrows),encoding='utf-8')
print(json.dumps(summary,indent=2))
