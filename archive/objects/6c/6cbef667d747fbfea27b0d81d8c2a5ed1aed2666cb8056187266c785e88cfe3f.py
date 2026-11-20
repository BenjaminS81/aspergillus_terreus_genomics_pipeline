"""User-authorised cold recovery of the paused genome, then sequential continuation."""
import sys,subprocess,json,fcntl,shutil
from pathlib import Path
import worker as w
def protected_current():
 units=json.loads(subprocess.check_output(['systemctl','--user','list-units','--type=service','--state=running','--output=json','--no-pager'],text=True))
 names=[u['unit'] for u in units if u['unit'].startswith('kleb-global-100gb-')]
 records=[]
 for name in names:
  p=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',name,'-p','MemoryMax','-p','MemoryCurrent','-p','ActiveState'],text=True).splitlines())
  assert p['MemoryMax']=='100000000000',('Protected limit changed',name,p)
  records.append({'unit':name,**p})
 return {'policy':'read-only discovery of current protected services; never alter them','services':records,'status':'verified' if records else 'no active Klebsiella 100 GB service found'}
w.protected=protected_current
def main():
 lock=(w.OUT/'execution.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 limits=w.limits();assert int(limits['memory.max'])==w.BUSCO_CAP
 pause=w.load(w.OUT/'PAUSED_BY_USER.json');folder=Path(pause['checkpoint']);snap=w.load(folder/'checkpoint.json')
 old=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',snap['unit'],'-p','MainPID','-p','FreezerState'],text=True).splitlines())
 assert old['MainPID']!=snap['service']['MainPID'],'Original paused worker still exists; use live resume'
 subprocess.run([sys.executable,str(w.OUT/'scripts/resume_checkpoint.py'),'--check'],check=True)
 assert w.load(w.OUT/'reports/annotation_acceptance.json')['sample_labels']==77
 for name in ['annotation_boundary_resolution','annotation_CDS_validation','orthofinder_input_recovery','functional_reanalysis']:
  assert w.load(w.OUT/'reports'/(name+'.json'))['exit_code']==0,('Previous stage not complete',name)
 job=snap['current_job'];queue=w.load(w.OUT/'data/genome_busco_execution_queue.json');item=next(j for j in queue if j['job']==job)
 argv=snap['cold_restart_argv'];assert '--restart' in argv and '--force' not in argv
 inp=w.lp(item['input_windows']);assert w.sha(inp)==item['expected_file_sha256']
 lineage=w.load(w.PREP/'data/busco_lineage_lock.json');base=w.lp(lineage['path'])
 print('Verifying lineage and retained checkpoint before cold resume',flush=True)
 w.save(w.OUT/'status.json',{'state':'preflight','phase':'cold_resume','current_job':job,'at':w.now()})
 for x in lineage['files']:assert w.sha(base/x['path'].replace('\\','/'))==x['sha256']
 archive=folder/'cold_recovery';archive.mkdir(exist_ok=False)
 for p in [w.OUT/'logs'/(job+'.log'),w.OUT/'reports'/(job+'.json')]:shutil.copy2(p,archive/p.name)
 (w.OUT/'PAUSED_BY_USER.json').rename(archive/'PAUSE_RESOLVED.json')
 w.save(archive/'resume_authorization.json',{'at':w.now(),'instruction':'Resume the Job where you paused','mode':'disk checkpoint recovery','completed_genomes_retained':snap['completed_genomes'],'limits':limits,'protected':protected_current(),'argv':argv})
 prefix=Path(argv[argv.index('--out_path')+1])/argv[argv.index('--out')+1]
 print('RESUME',job,'using --restart; completed stages and 28 validated genomes retained',flush=True)
 result=w.run(job,argv,item['environment']);assert result['exit_code']==0,('BUSCO recovery failed',result['exit_code'])
 w.save(w.OUT/'reports'/(job+'.validated.json'),w.check_busco(prefix,inp,item['expected_file_sha256']))
 print('Recovered genome validated; continuing remaining queue',flush=True)
 w.busco();w.protein_busco();w.core()
 for name,exe,script in [('supported_core_RF_NMDS','/home/benjamin/miniforge3/envs/klebsiella_analysis/bin/python','compare_supported_core.py'),('read_QC_and_mapping',sys.executable,'read_qc.py'),('gate_acceptance_report',sys.executable,'acceptance_report.py')]:
  r=w.run(name,[exe,str(w.OUT/'scripts'/script)]);assert r['exit_code']==0,(name,'failed; inspect log')
 w.save(w.OUT/'status.json',{'state':'stage_complete','phase':'cold_resume_and_continuation','finished':w.now(),'gate_closure':'requires integrated scientific acceptance'})
if __name__=='__main__':
 try:main()
 except Exception as e:
  w.save(w.OUT/'status.json',{'state':'stopped','phase':'cold_resume','at':w.now(),'reason':str(e)});raise
