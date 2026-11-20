"""Sequential, resource-limited execution with per-input provenance and result gates."""
from pathlib import Path
import argparse,os,sys,json,hashlib,subprocess,time,datetime,signal,fcntl,shutil
OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
PREP=ROOT/'gate_preparation/2026-11-18_v1'
PROTECTED='kleb-global-100gb-joint7830recovery1-20261118t003000z.service'
CAP=4*1024**3;BUSCO_CAP=20_000_000_000//4096*4096;RESERVE=12*1024**3
for name in ['logs','reports','data','runs']:(OUT/name).mkdir(parents=True,exist_ok=True)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,v):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.pending');tmp.write_text(json.dumps(v,indent=2));tmp.replace(p)
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def lp(p):
 s=str(p).replace('\\','/')
 return Path('/mnt/'+s[0].lower()+s[2:]) if len(s)>2 and s[1]==':' else Path(s)
def available():return int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
def protected():
 r=subprocess.run(['systemctl','--user','show',PROTECTED,'-p','MemoryMax','-p','MemoryCurrent','-p','ActiveState','-p','ControlGroup'],capture_output=True,text=True,check=True)
 d=dict(x.split('=',1) for x in r.stdout.splitlines() if '=' in x)
 if d.get('MemoryMax')!='100000000000':raise RuntimeError('Protected service limit differs from recorded 100 GB; no job launched')
 return d
def limits():
 cg=Path('/sys/fs/cgroup')/next(x.split(':',2)[2] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::')).lstrip('/')
 if not cg.name.startswith('terreus-gates-'):raise RuntimeError('Wrong cgroup')
 d={k:(cg/k).read_text().strip() for k in ['memory.max','memory.high','memory.swap.max','cpu.max']}
 q,p=d['cpu.max'].split()
 if d['memory.max']=='max' or int(d['memory.max'])>BUSCO_CAP or d['memory.swap.max']!='0' or q=='max' or int(q)>int(p):raise RuntimeError('Hard limits not applied')
 return {'cgroup':str(cg),**d}
def set_stage_limit(job):
 current=limits();unit=Path(current['cgroup']).name
 cap=BUSCO_CAP if job.startswith('BUSCO_') else CAP
 high=18_000_000_000//4096*4096 if job.startswith('BUSCO_') else 3*1024**3
 assert unit.startswith('terreus-gates-') and unit.endswith('.service')
 subprocess.run(['systemctl','--user','set-property','--runtime',unit,'MemoryMax='+str(cap),'MemoryHigh='+str(high),'MemorySwapMax=0'],check=True)
 applied=limits();assert int(applied['memory.max'])==cap and int(applied['memory.high'])==high
 return cap

def run(job,argv,env=None):
 cap=set_stage_limit(job)
 while available()<RESERVE+cap:
  save(OUT/'status.json',{'state':'waiting_for_memory_headroom','next_job':job,'at':now(),'available_bytes':available(),'required_available_bytes':RESERVE+cap,'active_analysis_jobs':0})
  time.sleep(30)
 record={'job':job,'started':now(),'argv':argv,'limits':limits(),'protected_before':protected(),'available_before':available()}
 log=OUT/'logs'/(job+'.log');status=OUT/'reports'/(job+'.json')
 save(OUT/'status.json',{'state':'running','current_job':job,'started':record['started'],'resource_policy':record['limits']})
 e=os.environ.copy();e.update(env or {})
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:e[k]='1'
 with log.open('w') as out:
  child=subprocess.Popen(argv,env=e,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
  record['pid']=child.pid;save(status,record)
  while child.poll() is None:
   if available()<RESERVE:
    os.killpg(child.pid,signal.SIGTERM)
    try:child.wait(timeout=5)
    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
    record['resource_guard']='Stopped this owned job because available RAM fell below 12 GiB'
    break
   time.sleep(1)
 record.update({'finished':now(),'exit_code':child.returncode,'available_after':available(),'protected_after':protected()})
 save(status,record)
 if 'resource_guard' in record:raise RuntimeError(record['resource_guard'])
 return record
def check_busco(outdir,input_path,digest,mode='genome'):
 ss=list(outdir.glob('short_summary.specific.*.json'))
 assert len(ss)==1,('summary count',len(ss))
 d=load(ss[0]);r=d['results'];pa=d['parameters'];lin=d['lineage_dataset']
 assert lin['name']=='eurotiales_odb12.2' and int(r['n_markers'])==4240
 assert (('genome' in pa['mode']) if mode=='genome' else (pa['mode']=='proteins')) and d['versions']['busco']=='6.1.0'
 assert r['Complete BUSCOs']==r['Single copy BUSCOs']+r['Multi copy BUSCOs']
 assert r['Complete BUSCOs']+r['Fragmented BUSCOs']+r['Missing BUSCOs']==4240
 table=outdir/'run_eurotiales_odb12.2/full_table.tsv'
 ids={line.split('\t',1)[0] for line in table.read_text().splitlines() if line and not line.startswith('#')}
 assert len(ids)==4240,('BUSCO full table marker coverage',len(ids))
 assert sha(input_path)==digest,'Input changed during run'
 return {'summary':str(ss[0]),'summary_sha256':sha(ss[0]),'full_table_sha256':sha(table),'full_table_unique_markers':len(ids),'input_sha256':digest,'results':r,'versions':d['versions'],'lineage':lin,'validation':'PASS'}
def busco(pilot=False):
 queue=load(PREP/'data/busco_missing_genome_queue.json')
 metrics=load(ROOT/'comparative_visualization/2026-11-13/data/genome_metrics.json')
 queue.sort(key=lambda j:metrics[j['job'].removeprefix('BUSCO_genome_')]['contigs'])
 if pilot:queue=queue[:1]
 else:
  # Old summaries lack run-time input hashes. Revalidate the remaining unique
  # assemblies after the 35 missing cases, so the final panel is fully pinned.
  import copy
  existing={j['job'] for j in queue};template=queue[0]
  for a in load(PREP/'data/canonical_assemblies.json'):
   name='BUSCO_genome_'+a['representative']
   if name in existing:continue
   j=copy.deepcopy(template);j.update({'job':name,'input_windows':a['windows_path'],'expected_file_sha256':a['current_file_sha256'],'reason':'historical input hash unavailable; prospective provenance verification'})
   j['argv'][j['argv'].index('--in')+1]=a['subject'];j['argv'][j['argv'].index('--out')+1]='genome_'+a['representative'];queue.append(j)
  assert len(queue)==72
  save(OUT/'data/genome_busco_execution_queue.json',queue)
 lock=load(PREP/'data/busco_lineage_lock.json');lineage=lp(lock['path'])
 print('Verifying pinned BUSCO lineage',flush=True)
 for item in lock['files']:assert sha(lineage/item['path'].replace('\\','/'))==item['sha256'],item['path']
 save(OUT/'reports/lineage_verification.json',{'completed':now(),'files':len(lock['files']),'lock_sha256':sha(PREP/'data/busco_lineage_lock.json'),'status':'PASS'})
 for i,j in enumerate(queue,1):
  job=j['job'];inp=lp(j['input_windows']);digest=j['expected_file_sha256'];argv=j['argv'][:]
  assert sha(inp)==digest,'Input changed: '+job
  argv[argv.index('--out_path')+1]=str(OUT/'runs/busco')
  outdir=OUT/'runs/busco'/argv[argv.index('--out')+1]
  done=OUT/'reports'/(job+'.validated.json')
  if done.exists():
   verified=check_busco(outdir,inp,digest);save(done,verified);print('Verified completed',job,flush=True);continue
  if outdir.exists():raise RuntimeError('Incomplete output requires inspection, not automatic overwrite: '+str(outdir))
  (OUT/'runs/busco').mkdir(exist_ok=True)
  print(f'RUN {i}/{len(queue)} {job}',flush=True)
  rec=run(job,argv,j['environment'])
  if rec['exit_code']!=0:raise RuntimeError('BUSCO failed: '+job)
  verified=check_busco(outdir,inp,digest);save(done,verified)
  print('PASS',job,verified['results'].get('one_line_summary'),flush=True)
 save(OUT/'reports/busco_stage.json',{'finished':now(),'pilot':pilot,'validated':len(list((OUT/'reports').glob('BUSCO_genome_*.validated.json'))),'planned_unique_jobs':1 if pilot else 72,'original_missing_jobs':35})
def protein_busco():
 rows=load(OUT/'reports/annotation_release_summary.json');groups={}
 for row in rows:
  if row.get('primary_genes',0)==0:continue
  p=OUT/'runs/annotation_release'/row['sample']/'primary_validated.faa';digest=sha(p)
  groups.setdefault(digest,{'input':str(p),'samples':[],'source_model_statuses':[]})['samples'].append(row['sample'])
  groups[digest]['source_model_statuses'].append(row['technical_status'])
 save(OUT/'data/protein_busco_groups.json',groups)
 env='/home/benjamin/.local/share/Achromobacter_ANI/envs/busco';lineage=lp(load(PREP/'data/busco_lineage_lock.json')['path'])
 for digest,g in groups.items():
  sample=g['samples'][0];job='BUSCO_proteins_'+sample;outparent=OUT/'runs/busco_proteins';outparent.mkdir(exist_ok=True);outdir=outparent/('proteins_'+sample);done=OUT/'reports'/(job+'.validated.json')
  if done.exists():check_busco(outdir,g['input'],digest,'proteins');continue
  if outdir.exists():raise RuntimeError('Incomplete protein BUSCO requires review: '+str(outdir))
  argv=[env+'/bin/python3.12',env+'/bin/busco','--in',g['input'],'--out','proteins_'+sample,'--out_path',str(outparent),'--mode','proteins','--lineage_dataset',str(lineage),'--offline','--cpu','1','--opt-out-run-stats']
  rec=run(job,argv,{'PATH':env+'/bin:/usr/local/bin:/usr/bin:/bin'})
  if rec['exit_code']!=0:raise RuntimeError('Protein BUSCO failed: '+sample)
  result=check_busco(outdir,g['input'],digest,'proteins');result.update(g);result['scope']='validated primary gene set; excludes rejected CDS models and is separate from genome mode';save(done,result)
  print('PASS',job,result['results'].get('one_line_summary'),flush=True)
 save(OUT/'reports/protein_busco_stage.json',{'finished':now(),'unique_gene_sets':len(groups),'sample_labels':sum(len(g['samples']) for g in groups.values())})
def annotations():
 r=run('annotation_boundary_resolution',[sys.executable,str(OUT/'scripts/resolve_annotation_links.py')])
 if r['exit_code']!=0:raise RuntimeError('Annotation boundary resolver failed; inspect log')
 r=run('annotation_CDS_validation',[sys.executable,str(OUT/'scripts/validate_annotations.py')])
 if r['exit_code']!=0:raise RuntimeError('Annotation validator failed; inspect log')
def core():
 plan=load(PREP/'data/core_support_job_plan.json');inp=lp(plan['input_windows'])
 assert sha(inp)==plan['input_sha256']
 exe='/home/benjamin/miniforge3/envs/klebsiella_analysis/bin/iqtree3'
 version=subprocess.run([exe,'--version'],capture_output=True,text=True,check=True).stdout
 assert 'version 3.1.1 ' in version
 original=load(ROOT/'review_reports/2026-11-10/audit_evidence.json')['trees']['genome']
 assert sha(Path(str(inp)+'.treefile'))==original['sha256'],'Starting core tree changed'
 prefix=OUT/'runs/core_support/core_GTR_F_R10';prefix.parent.mkdir(exist_ok=True)
 argv=[exe,'-s',str(inp),'-m','GTR+F+R10','-t',str(inp)+'.treefile','-B','1000','-alrt','1000','-bnni','-wbt','-seed','20261118','-nt','1','-mem','2G','--prefix',str(prefix)]
 if Path(str(prefix)+'.iqtree').exists() and (OUT/'reports/core_support.validated.json').exists():return
 r=run('core_support',argv)
 if r['exit_code']!=0:
  save(OUT/'reports/core_support.blocked.json',{'exit_code':r['exit_code'],'reason':'Inspect core_support.log; no automatic cap increase','resource_cap':CAP});return
 from Bio import Phylo
 tree=Phylo.read(str(prefix)+'.treefile','newick')
 names={t.name for t in tree.get_terminals()}
 expected=set(load(ROOT/'review_reports/2026-11-10/audit_evidence.json')['trees']['genome']['tips'])
 assert names==expected
 boot=Path(str(prefix)+'.ufboot');n=sum(line.strip().endswith(';') for line in boot.open())
 assert n>=1000
 supports=[x.name for x in tree.get_nonterminals() if x.name and '/' in x.name]
 assert supports,'SH-aLRT/UFBoot support labels absent'
 save(OUT/'reports/core_support.validated.json',{'status':'PASS','alignment_sha256':plan['input_sha256'],'version':version.strip(),'tree_sha256':sha(Path(str(prefix)+'.treefile')),'tips':len(names),'bootstrap_trees':n,'supported_internal_nodes':len(supports),'rooting_and_model_sensitivity':'still required; support alone does not close G3'})
def main():
 p=argparse.ArgumentParser();p.add_argument('phase',choices=['pilot','busco','annotations','full','resume_all']);a=p.parse_args()
 lock=(OUT/'execution.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 save(OUT/'reports/resource_preflight.json',{'at':now(),'limits':limits(),'protected':protected(),'memory_available':available()})
 try:
  save(OUT/'status.json',{'state':'preflight','phase':a.phase,'started':now()})
  if a.phase=='resume_all':busco(pilot=True)
  if a.phase=='annotations':annotations()
  elif a.phase in {'full','resume_all'}:
   annotations()
   r=run('orthofinder_input_recovery',[sys.executable,str(OUT/'scripts/recover_orthofinder_inputs.py')])
   if r['exit_code']!=0:raise RuntimeError('Orthofinder input recovery failed')
   r=run('functional_reanalysis',[sys.executable,str(OUT/'scripts/reconcile_functional.py')])
   if r['exit_code']!=0:raise RuntimeError('Functional parser failed; inspect log')
   subprocess.run([sys.executable,str(OUT/'scripts/acceptance_report.py')],check=True)
   busco()
   protein_busco()
   core()
   r=run('supported_core_RF_NMDS',['/home/benjamin/miniforge3/envs/klebsiella_analysis/bin/python',str(OUT/'scripts/compare_supported_core.py')])
   if r['exit_code']!=0:raise RuntimeError('Updated core comparison failed; inspect log')
   r=run('read_QC_and_mapping',[sys.executable,str(OUT/'scripts/read_qc.py')])
   if r['exit_code']!=0:raise RuntimeError('Read QC stage failed; inspect log')
   r=run('gate_acceptance_report',[sys.executable,str(OUT/'scripts/acceptance_report.py')])
   if r['exit_code']!=0:raise RuntimeError('Acceptance report failed')
  else:busco(a.phase=='pilot')
  save(OUT/'status.json',{'state':'stage_complete','phase':a.phase,'finished':now(),'gate_closure':'requires integrated acceptance; not implied by process completion'})
 except Exception as e:
  save(OUT/'status.json',{'state':'stopped','phase':a.phase,'at':now(),'reason':str(e)})
  reporter=OUT/'scripts/acceptance_report.py'
  if reporter.exists():subprocess.run([sys.executable,str(reporter)],check=False)
  raise
if __name__=='__main__':main()
