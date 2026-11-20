"""Sequential, streaming retrieval of two approved genome records; no analysis jobs."""
from pathlib import Path
import json,urllib.request,hashlib,gzip,datetime,time,os,ctypes
ROOT=Path(__file__).resolve().parent
AUDIT=ROOT/'review_reports/2026-11-20_public_genome_audit_v1'
OUT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\data\external\public_additions\2026-11-20_v1')
REVIEW=ROOT/'review_reports/2026-11-20_public_genome_audit_v3_fetch';REVIEW.mkdir(exist_ok=True)
OUT.mkdir(parents=True,exist_ok=True)
def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
if os.name=='nt':
 k=ctypes.windll.kernel32;k.GetCurrentProcess.restype=ctypes.c_void_p
 handle=k.GetCurrentProcess();k.SetPriorityClass.argtypes=[ctypes.c_void_p,ctypes.c_ulong];k.SetPriorityClass(handle,0x4000)
 proc=ctypes.c_size_t();system=ctypes.c_size_t();k.GetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.POINTER(ctypes.c_size_t),ctypes.POINTER(ctypes.c_size_t)]
 if k.GetProcessAffinityMask(handle,ctypes.byref(proc),ctypes.byref(system)):
  k.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t];k.SetProcessAffinityMask(handle,proc.value & -proc.value)
manifest={'started_utc':now(),'destination':str(OUT),'policy':'Sequential single-thread streaming, 1 MiB buffers; below-normal priority and one CPU affinity on Windows. No changes to existing analyses or queues.','assemblies':[],'existing_read_aliases':[],'status':'running'}
def checkpoint():
 save(OUT/'download_manifest.json',manifest);save(REVIEW/'download_manifest.json',manifest)
def fetch(url,dest):
 if dest.exists():return
 partial=dest.with_name(dest.name+'.part')
 for attempt in range(3):
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Terreus-genome-review/1.0'}),timeout=60) as r,partial.open('wb') as f:
    while block:=r.read(1024*1024):f.write(block)
   partial.replace(dest);return
  except Exception:
   if attempt==2:raise
   time.sleep(2)
def hashes(p):
 md5=hashlib.md5();sha=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024*1024):md5.update(b);sha.update(b)
 return {'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'bytes':p.stat().st_size}
reports=json.loads((AUDIT/'ncbi_assemblies.json').read_text());summ=json.loads((AUDIT/'ncbi_assembly_summaries.json').read_text())['result']
targets=['GCA_061584345.1','GCA_988135485.1']
try:
 for accession in targets:
  record=next(x for x in reports if x['accession']==accession);summary=next(summ[u] for u in summ['uids'] if summ[u]['assemblyaccession']==accession)
  dest=OUT/'assemblies'/accession;dest.mkdir(parents=True,exist_ok=True)
  base=summary['ftppath_genbank'].replace('ftp://','https://');stem=base.rsplit('/',1)[-1]
  item={'accession':accession,'source_directory':base,'biosample':record['assembly_info']['biosample']['accession'],'role':'new genome candidate' if accession==targets[0] else 'provenance reconciliation with existing ERR17672364/T59-F representations','files':[],'status':'downloading'}
  manifest['assemblies'].append(item);checkpoint();print('Downloading',accession,flush=True)
  save(dest/'ncbi_dataset_report.json',record)
  fetch(base+'/md5checksums.txt',dest/'md5checksums.txt')
  published={}
  for line in (dest/'md5checksums.txt').read_text().splitlines():
   parts=line.split(maxsplit=1)
   if len(parts)==2:published[parts[1].lstrip('./')]=parts[0]
  names=[stem+'_genomic.fna.gz',stem+'_assembly_report.txt',stem+'_assembly_stats.txt']
  for name in names:
   if name not in published:
    if name.endswith('_genomic.fna.gz'):raise RuntimeError('Required genomic FASTA missing from published checksum list')
    continue
   path=dest/name;url=base+'/'+name;fetch(url,path);h=hashes(path)
   if h['md5']!=published[name]:raise RuntimeError('Published MD5 mismatch: '+str(path))
   item['files'].append({'path':str(path),'url':url,'published_md5':published[name],'verified_utc':now(),**h})
  packed=dest/(stem+'_genomic.fna.gz');fasta=dest/(stem+'_genomic.fna')
  part=fasta.with_name(fasta.name+'.part')
  with gzip.open(packed,'rb') as z,part.open('wb') as f:
   while block:=z.read(1024*1024):f.write(block)
  part.replace(fasta)
  lengths=[];length=None;ns=0
  with fasta.open('rb') as f:
   for line in f:
    if line.startswith(b'>'):
     if length is not None:lengths.append(length)
     length=0
    else:
     if length is None:raise RuntimeError('Invalid FASTA header')
     b=line.strip();length+=len(b);ns+=b.upper().count(b'N')
  if length is not None:lengths.append(length)
  expected=int(record['assembly_stats']['total_sequence_length']);assert sum(lengths)==expected,(accession,sum(lengths),expected)
  item['fasta']={'path':str(fasta),**hashes(fasta),'sequences':len(lengths),'bases':sum(lengths),'N_bases':ns,'total_length_matches_ncbi':True,'gzip_integrity_passed':True}
  item['status']='downloaded_and_verified';checkpoint();print('Verified',accession,len(lengths),'sequences',sum(lengths),'bases',flush=True)
 # Existing anonymous FASTQ aliases: check integrity and record/base counts without downloading copies.
 ena=json.loads((AUDIT/'ena_read_run.json').read_text())
 pairs=[('SRR30681804',Path(r'D:\Illumina\David_Eisele\terreus_genotyping\inbox\ATSec\fastq\ATsec\ATsec.fastq.gz')),('SRR30681805',Path(r'D:\Illumina\David_Eisele\terreus_genotyping\inbox\ATSec\fastq\wt\wt.fastq.gz'))]
 for run,path in pairs:
  r=next(x for x in ena if x['run_accession']==run)
  item={'public_run':run,'existing_path':str(path),'sample_title':r['sample_title'],'public_read_count':int(r['read_count']),'public_base_count':int(r['base_count']),'status':'validating_existing_alias','exact_public_sequence_identity':'not established by counts alone'}
  manifest['existing_read_aliases'].append(item);checkpoint();print('Validating existing read alias',run,str(path),flush=True)
  count=0;bases=0;first_headers=[]
  with gzip.open(path,'rb') as f:
   while True:
    header=f.readline()
    if not header:break
    seq=f.readline().rstrip(b'\r\n');plus=f.readline();qual=f.readline().rstrip(b'\r\n')
    if not header.startswith(b'@') or not plus.startswith(b'+') or len(seq)!=len(qual):raise RuntimeError('FASTQ structural error at '+str(count))
    count+=1;bases+=len(seq)
    if len(first_headers)<3:first_headers.append(header.decode(errors='replace').strip())
    if count%100000==0:print(run,count,'reads checked',flush=True)
  item.update({'local_read_count':count,'local_base_count':bases,'compressed_bytes':path.stat().st_size,'gzip_and_fastq_structure_valid':True,'first_headers':first_headers,'counts_match_public':count==int(r['read_count']) and bases==int(r['base_count']),'status':'existing_data_validated','download_decision':'Existing WT/sector data retained; no duplicate public read download. Accession alias remains provisional unless exact sequence identity is established.'})
  checkpoint();print('Existing read alias result',run,count,bases,item['counts_match_public'],flush=True)
 manifest['status']='complete';manifest['finished_utc']=now();checkpoint()
 text='Relevant public genome additions — 20 November 2026\n\nDownloaded and verified: CMI 44339 (GCA_061584345.1), 59TF (GCA_988135485.1).\nGenomic FASTA gzip and extracted FASTA, assembly reports/statistics, NCBI metadata and published MD5 checksums are retained. Each extracted FASTA total length matches NCBI.\n\nThese are staged additions, not yet included in any tree or active QC queue. CMI 44339 needs standard genome/locus QC; 59TF needs provenance reconciliation with ERR17672364 and T59-F before independent inclusion.\n\nExisting WT/sector FASTQs were found in inbox/ATSec and validated by streaming gzip/FASTQ and read/base counts. Public accession identity is not proven solely by matching counts. No duplicate read download was made. See download_manifest.json for exact findings and paths.\n\nExcluded: DAH1047FM (347 bp, not a usable genome), historical ATCC20542 partial record (138.5 kb), legacy capillary runs. Already downloaded T3_Kankrej/HNGD-TM15 were not downloaded again. Existing analyses and resource caps were unchanged.\n'
 (OUT/'README.txt').write_text(text,encoding='utf-8');(REVIEW/'README.txt').write_text(text,encoding='utf-8')
 print('COMPLETE',str(OUT),flush=True)
except Exception as ex:
 manifest['status']='failed';manifest['error']=repr(ex);checkpoint();raise
