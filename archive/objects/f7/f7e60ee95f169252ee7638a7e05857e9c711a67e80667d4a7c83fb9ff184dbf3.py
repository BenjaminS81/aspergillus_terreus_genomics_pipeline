"""Sequential raw-pair QC and assembly mapping; source reads/assemblies stay read-only."""
from pathlib import Path
import re,gzip,json,collections,hashlib,subprocess,shutil,os
from worker import OUT,PREP,lp,load,save,sha,now
FASTP=OUT/'tools/fastp.1.4.0';ENV='/home/benjamin/miniforge3/envs/klebsiella_analysis/bin'
FASTP_SHA='e56c3218b76cc1d6d8a7a32315fa6a41407f5d271fe561fc62cee15941190844'
MATE=re.compile(r'^(.*?)([_.]R?)([12])((?:_\d+)?\.(?:fastq|fq)(?:\.gz)?)$')
def command(argv,log,stdout=None):
 with Path(log).open('w') as err:
  if stdout:
   with Path(stdout).open('w') as out:r=subprocess.run(list(map(str,argv)),stdout=out,stderr=err)
  else:r=subprocess.run(list(map(str,argv)),stdout=err,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError('Exit '+str(r.returncode)+': '+str(argv[0]))
def first_headers(p):
 op=gzip.open if str(p).endswith('.gz') else open
 result=[]
 with op(p,'rt') as f:
  for i in range(1000):
   h=f.readline()
   if not h:break
   seq=f.readline().strip();plus=f.readline();qual=f.readline().strip()
   if not h.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):raise ValueError('Malformed FASTQ in header precheck')
   result.append(re.sub(r'/[12]$','',h[1:].split()[0]))
 return result
def main():
 assert sha(FASTP)==FASTP_SHA
 manifest=load(PREP/'data/analysis_eligibility.json');groups={};unresolved=[]
 for m in manifest:
  pairs=collections.defaultdict(dict)
  for raw in m['raw_read_candidates']:
   p=lp(raw);hit=MATE.match(p.name)
   if not hit:continue
   key=(str(p.parent),hit[1],hit[2],hit[4]);pairs[key][hit[3]]=str(p)
  for key,pair in pairs.items():
   if set(pair)!={'1','2'}:unresolved.append({'sample':m['sample'],'reason':'incomplete_mate_pair','files':pair});continue
   k=(pair['1'],pair['2']);groups.setdefault(k,[]).append(m)
  if m['raw_read_candidates'] and not pairs:unresolved.append({'sample':m['sample'],'reason':'unrecognized_mate_filenames'})
 save(OUT/'data/read_pair_groups.json',[{'read1':k[0],'read2':k[1],'candidate_sample_labels':[m['sample'] for m in ms]} for k,ms in groups.items()])
 save(OUT/'reports/read_pair_exceptions.json',unresolved)
 versions={}
 for name,exe in [('fastp',str(FASTP)),('minimap2',ENV+'/minimap2'),('samtools',ENV+'/samtools')]:
  r=subprocess.run([exe,'--version'],capture_output=True,text=True,check=True);versions[name]=(r.stdout+r.stderr).strip()
 save(OUT/'data/read_qc_tool_versions.json',versions)
 results=[]
 for pair,samples in groups.items():
  key=hashlib.sha256('\n'.join(pair).encode()).hexdigest()[:16];folder=OUT/'runs/read_qc'/key;folder.mkdir(parents=True,exist_ok=True)
  report=folder/'validation.json';save(OUT/'reports/read_qc_progress.json',{'state':'pair_QC','pair':key,'sample_candidates':[m['sample'] for m in samples],'at':now()})
  try:
   if report.exists():
    rec=load(report)
    assert all(Path(p).stat().st_size==r['bytes'] and Path(p).stat().st_mtime_ns==r['mtime_ns'] for p,r in rec['inputs'].items()),'Source changed after QC'
   else:
    hh=[first_headers(p) for p in pair];assert hh[0] and hh[0]==hh[1],'Mates disagree in the first 1000 records'
    inputs={p:{'sha256':sha(p),'bytes':Path(p).stat().st_size,'mtime_ns':Path(p).stat().st_mtime_ns} for p in pair}
    args=[str(FASTP),'-i',pair[0],'-I',pair[1],'-w','1','--detect_adapter_for_pe','--dont_eval_duplication','--qualified_quality_phred','15','--unqualified_percent_limit','40','--n_base_limit','5','--length_required','15','--json',str(folder/'fastp.json'),'--html',str(folder/'fastp.html')]
    command(args,folder/'fastp.log');qc=load(folder/'fastp.json')
    assert qc['summary']['before_filtering']['total_reads']>0
    assert all(Path(p).stat().st_size==r['bytes'] and Path(p).stat().st_mtime_ns==r['mtime_ns'] for p,r in inputs.items())
    rec={'pair':key,'inputs':inputs,'candidate_sample_labels':[m['sample'] for m in samples],'header_pair_precheck_records':len(hh[0]),'argv':args,'summary':qc['summary'],'status':'PASS','policy':'No trimmed FASTQ files emitted. Before-filtering values describe raw reads; after-filtering values describe this QC policy, not historical assembly inputs. No duplicate estimate requested.'};save(report,rec)
   results.append(rec);save(OUT/'reports/read_qc_summary.json',results)
   mapped=set()
   for m in samples:
    assembly_hash=m['sequence_sha256']
    if assembly_hash in mapped:continue
    mapped.add(assembly_hash);dest=OUT/'runs/read_mapping'/(key+'_'+assembly_hash[:12]);dest.mkdir(parents=True,exist_ok=True);done=dest/'validation.json'
    if done.exists():continue
    read_size=sum(r['bytes'] for r in rec['inputs'].values())
    if shutil.disk_usage(OUT).free<3*read_size+64*1024**3:raise RuntimeError('Insufficient reserved output disk space for BAM')
    genome=lp(m['assembly']);assert sha(genome)==m['assembly_file_sha256'];bam=dest/'reads.sorted.bam'
    assert not bam.exists(),'Incomplete BAM needs inspection before retry'
    ma=[ENV+'/minimap2','-ax','sr','--secondary=no','-t','1','-I','1G',str(genome),*pair]
    sa=[ENV+'/samtools','sort','-@','0','-m','256M','-T',str(dest/'sort_tmp'),'-o',str(bam),'-']
    save(OUT/'reports/read_qc_progress.json',{'state':'mapping','pair':key,'assembly_sample':m['sample'],'at':now()})
    with (dest/'minimap2.log').open('w') as ml,(dest/'sort.log').open('w') as sl:
     first=subprocess.Popen(ma,stdout=subprocess.PIPE,stderr=ml)
     second=subprocess.Popen(sa,stdin=first.stdout,stdout=sl,stderr=subprocess.STDOUT);first.stdout.close()
     rc2=second.wait();rc1=first.wait()
     assert rc1==0 and rc2==0,('Mapping/sort failed',rc1,rc2)
    command([ENV+'/samtools','quickcheck','-v',bam],dest/'quickcheck.log')
    command([ENV+'/samtools','index','-@','0',bam],dest/'index.log')
    command([ENV+'/samtools','flagstat','-@','0','-O','json',bam],dest/'flagstat.log',dest/'flagstat.json')
    command([ENV+'/samtools','coverage','--ff','0xF04','-q','20','-Q','20','-o',dest/'coverage.tsv',bam],dest/'coverage.log')
    save(done,{'status':'PASS','read_pair':key,'read_input_hashes':rec['inputs'],'assembly_file_sha256':m['assembly_file_sha256'],'assembly_sequence_sha256':assembly_hash,'mapped_sample_labels':[x['sample'] for x in samples if x['sequence_sha256']==assembly_hash],'mapping_argv':ma,'sort_argv':sa,'bam_sha256':sha(bam),'coverage_filters':{'minimum_mapping_quality':20,'minimum_base_quality':20,'excluded_flags':'0xF04'},'duplicate_policy':'no duplicate marking/removal; depth is read-depth QC, not validated CNV','sample_identity':'candidate name/accession link; mapping is not proof of culture equivalence'})
   print('Read QC/mapping complete',key,flush=True)
  except Exception as e:
   result={'pair':key,'status':'FAILED_OR_UNRESOLVED','reason':str(e),'candidate_sample_labels':[m['sample'] for m in samples]};results.append(result);save(folder/'exception.json',result);print('Read case unresolved',key,str(e),flush=True)
  save(OUT/'reports/read_qc_summary.json',results)
 save(OUT/'reports/read_qc_stage.json',{'finished':now(),'paired_file_groups':len(groups),'mapping_results':len(list((OUT/'runs/read_mapping').glob('*/validation.json'))) if (OUT/'runs/read_mapping').exists() else 0,'unresolved_filename_groups':len(unresolved),'case_status_counts':dict(collections.Counter(r['status'] for r in results)),'historical_downsampling_and_assembly_log_provenance':'not reconstructed by rerunning read QC'})
if __name__=='__main__':main()
