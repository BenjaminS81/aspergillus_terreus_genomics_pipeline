"""Read-only genome inventory audit; saves public metadata and review files only."""
from pathlib import Path
import json,urllib.request,urllib.parse,time,hashlib,datetime,subprocess,re,collections
ROOT=Path(__file__).resolve().parent
PROJECT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
OUT=ROOT/'review_reports/2026-11-20_public_genome_audit_v1';RAW=OUT/'raw';RAW.mkdir(parents=True,exist_ok=True)
def save(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
manifest=json.loads((OUT/'retrieval_manifest.json').read_text()) if (OUT/'retrieval_manifest.json').exists() else []
def get(url,name):
 p=RAW/name
 if p.exists():return json.loads(p.read_text())
 for attempt in range(3):
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Terreus-public-inventory-audit/1.0'}),timeout=55) as r:data=r.read()
   obj=json.loads(data);p.write_bytes(data);manifest.append({'file':name,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()});save('retrieval_manifest.json',manifest);return obj
  except Exception:
   if attempt==2:raise
   time.sleep(2)
def fetch():
 records=[];page=0;token=''
 while True:
  query={'page_size':1000}
  if token:query['page_token']=token
  d=get('https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/33178/dataset_report?'+urllib.parse.urlencode(query),'ncbi_terreus_page_'+str(page)+'.json');records+=d.get('reports',[]);token=d.get('next_page_token','');page+=1
  if not token:assert len(records)==int(d['total_count']);break
 save('ncbi_assemblies.json',records);print('NCBI assemblies',len(records),flush=True)
 for result in ['assembly','read_run','analysis']:
  fields=get('https://www.ebi.ac.uk/ena/portal/api/returnFields?result='+result+'&format=json','ena_'+result+'_fields.json');available={r['columnId'] for r in fields}
  desired=['accession','assembly_accession','analysis_accession','study_accession','secondary_study_accession','sample_accession','secondary_sample_accession','run_accession','experiment_accession','scientific_name','tax_id','strain','isolate','assembly_name','assembly_level','version','wgs_set','first_public','last_updated','description','study_title','sample_title','analysis_type','library_strategy','library_source','library_layout','instrument_platform','instrument_model','read_count','base_count','fastq_ftp','fastq_md5','fastq_bytes','submitted_ftp','submitted_md5','submitted_bytes']
  q='tax_tree(33178)'
  if result=='read_run':q+=' AND library_strategy="WGS"'
  if result=='analysis':q+=' AND analysis_type="SEQUENCE_ASSEMBLY"'
  params={'result':result,'query':q,'fields':','.join(x for x in desired if x in available),'format':'json','limit':0}
  data=get('https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode(params),'ena_'+result+'.json');save('ena_'+result+'.json',data);print('ENA',result,len(data),flush=True)
 # Independently enumerate NCBI Assembly IDs, including non-current records if returned.
 params={'db':'assembly','term':'txid33178[Organism:exp]','retmode':'json','retmax':10000}
 s=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urllib.parse.urlencode(params),'ncbi_assembly_search.json')['esearchresult'];assert int(s['count'])==len(s['idlist'])
 if s['idlist']:
  time.sleep(.4);d=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urllib.parse.urlencode({'db':'assembly','id':','.join(s['idlist']),'retmode':'json'}),'ncbi_assembly_summaries.json');save('ncbi_assembly_summaries.json',d)
 print('NCBI Assembly search',s['count'],flush=True)
def local():
 roots=[PROJECT/x for x in ['data/assemblies','data/raw','data/external','data/reference','inbox'] if (PROJECT/x).exists()]
 r=subprocess.run(['rg','--files','--hidden','--no-ignore',*[str(x) for x in roots]],capture_output=True,text=True,check=True)
 assets=json.loads((ROOT/'taxonomy_reference_build/2026-11-10/analysis/assembly_asset_inventory.json').read_text())
 paths=sorted(set(r.stdout.splitlines())|{x['windows_path'] for x in assets if Path(x['windows_path']).is_file()});assemblies=[];reads=[]
 for s in paths:
  p=Path(s);name=p.name
  if re.search(r'\.(fna|fa|fasta)(\.gz)?$',name,re.I) and re.search(r'GC[AF]_\d+\.\d+',name):
   st=p.stat();assemblies.append({'path':s,'accession':re.search(r'GC[AF]_\d+\.\d+',name)[0],'bytes':st.st_size,'mtime_ns':st.st_mtime_ns,'nonempty':st.st_size>0})
  if re.search(r'\.(fastq|fq)(\.gz)?$',name,re.I):
   acc=re.findall(r'[SED]RR\d+',name)
   if acc:
    st=p.stat();reads.append({'path':s,'run_accession':acc[0],'bytes':st.st_size,'mtime_ns':st.st_mtime_ns})
 save('local_candidate_paths.json',paths)
 save('local_assembly_files.json',assemblies);save('local_raw_read_files.json',reads);save('local_scan.json',{'root':str(PROJECT),'scanned_files':len(paths),'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'accession_named_assembly_files':len(assemblies),'public_run_named_read_files':len(reads),'scope':'Fresh scan of data/assemblies, data/raw, data/external, data/reference and inbox plus existing 620-asset assembly inventory with current existence checks. Filename/size availability, not a fresh FASTQ integrity audit. Uncatalogued files elsewhere are not excluded.'})
 print('Local scan',len(paths),'files;',len(assemblies),'assembly files;',len(reads),'read files',flush=True)
if __name__=='__main__':fetch();local()
