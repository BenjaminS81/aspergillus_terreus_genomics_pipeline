"""Save an audit/resume snapshot after a user-requested cgroup freeze. Never resumes."""
from pathlib import Path
import json,subprocess,datetime,time,hashlib,shutil
OUT=Path(__file__).resolve().parents[1]
UNIT='terreus-gates-all20gb-20261118.service'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
props=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',UNIT,'-p','FreezerState','-p','MainPID','-p','ControlGroup','-p','MemoryMax','-p','MemoryCurrent','-p','CPUQuotaPerSecUSec'],text=True).splitlines())
assert props['FreezerState']=='frozen'
cg=Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
assert 'frozen 1' in (cg/'cgroup.events').read_text()
def usage():return int(dict(l.split() for l in (cg/'cpu.stat').read_text().splitlines())['usage_usec'])
first=usage();time.sleep(1);last=usage();assert first==last,'Frozen service still using CPU'
at=datetime.datetime.now(datetime.timezone.utc);dest=OUT/'checkpoints'/at.strftime('%Y%m%dT%H%M%SZ_user_pause');dest.mkdir(parents=True)
state=read(OUT/'status.json');job=state['current_job'];job_record=read(OUT/'reports'/(job+'.json'))
argv=job_record['argv'];run=Path(argv[argv.index('--out_path')+1])/argv[argv.index('--out')+1]
validated=[]
for p in sorted((OUT/'reports').glob('BUSCO_genome_*.validated.json')):
 validated.append({'job':p.name.removesuffix('.validated.json'),'validation_file':str(p),'sha256':sha(p),'input_sha256':read(p)['input_sha256']})
queue=read(OUT/'data/genome_busco_execution_queue.json');done={r['job'] for r in validated}
checkpoints=[]
for p in run.rglob('*'):
 if p.is_file() and ('checkpoint' in p.name or p.suffix=='.ckp'):
  target=dest/'busco_checkpoints'/p.relative_to(run);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
  checkpoints.append({'original':str(p),'saved':str(target),'sha256':sha(p)})
for src in [OUT/'status.json',OUT/'reports'/(job+'.json'),OUT/'logs'/(job+'.log'),OUT/'data/genome_busco_execution_queue.json',OUT/'reports/lineage_verification.json']:
 target=dest/src.relative_to(OUT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,target)
shutil.copytree(OUT/'scripts',dest/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
snapshot={'paused_at':at.isoformat(),'scope':'all live Terreus analysis processes found in the owned service; unrelated workloads untouched','unit':UNIT,'service':props,'frozen_verified':True,'cpu_usage_unchanged_over_one_second':True,'current_job':job,'current_stage':'miniprot alignment in progress','completed_genomes':len(validated),'planned_genomes':len(queue),'validated_results':validated,'remaining_genome_jobs':[r['job'] for r in queue if r['job'] not in done],'busco_disk_checkpoints':checkpoints,'cold_restart_argv':argv+['--restart'],'live_resume':'Verify resources and snapshot; thaw the same frozen service. This preserves the in-progress alignment.','cold_resume':'After WSL shutdown/reboot, live memory is lost. Verify inputs and retained output; BUSCO --restart can reuse completed stages, but unfinished miniprot alignment must rerun. The queue worker currently refuses incomplete outputs, so cold recovery must explicitly validate this job before continuing.','automatic_resume':False}
(dest/'checkpoint.json').write_text(json.dumps(snapshot,indent=2))
(OUT/'PAUSED_BY_USER.json').write_text(json.dumps({'at':at.isoformat(),'checkpoint':str(dest),'unit':UNIT,'resume_only_on_user_instruction':True},indent=2))
(OUT/'status.json').write_text(json.dumps({'state':'paused_by_user','paused_at':at.isoformat(),'current_job':job,'completed_genome_runs':len(validated),'total_genome_runs':len(queue),'freezer_state':'frozen','checkpoint':str(dest),'live_memory_retained_bytes':int(props['MemoryCurrent']),'automatic_resume':False},indent=2))
(dest/'RESUME_NOTES.txt').write_text('USER-REQUESTED PAUSE\n\nAll project analysis processes are frozen. No jobs will be launched by this snapshot.\n28 completed genome results and their validation records remain on disk.\nThe current BUSCO alignment is retained in live process memory.\n\nLIVE RESUME (only after user instruction):\nCheck current competing workloads and memory headroom, verify the same service\nis frozen and its MainPID matches checkpoint.json, clear PAUSED_BY_USER.json,\nrestore the saved running status, then run:\nsystemctl --user thaw '+UNIT+'\n\nAFTER WSL SHUTDOWN / REBOOT:\nThe live alignment progress will be lost. Completed results and on-disk BUSCO\ncheckpoints remain. Use the recorded cold_restart_argv with unchanged input\nhashes and a capped service; validate the current genome before continuing\nthe queue. Never force-overwrite its output. Do not launch a duplicate worker\nwhile the frozen service exists.\n')
(OUT/'STATUS.txt').write_text('PAUSED BY USER\n\nCompleted genome BUSCO runs: '+str(len(validated))+'/'+str(len(queue))+'\nCurrent frozen job: '+job+'\nCheckpoint: '+str(dest)+'\nNo automatic resume. Current alignment remains in RAM; a reboot requires restarting its unfinished stage.\n')
print(json.dumps({'checkpoint':str(dest),'completed':len(validated),'current':job,'frozen':True,'cpu_usage_stable':True,'disk_checkpoints':len(checkpoints),'memory_retained_bytes':int(props['MemoryCurrent'])},indent=2))
