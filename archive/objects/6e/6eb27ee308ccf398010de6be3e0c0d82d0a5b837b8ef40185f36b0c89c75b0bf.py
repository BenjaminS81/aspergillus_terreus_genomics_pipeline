"""Default dry run. Execute one prepared BUSCO job only inside verified limits."""
from pathlib import Path
import argparse, os, sys, json, hashlib, subprocess, datetime
OUT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser()
p.add_argument('job',help='exact job ID in busco_missing_genome_queue.json')
p.add_argument('--execute',action='store_true')
a=p.parse_args()
jobs=json.loads((OUT/'data/busco_missing_genome_queue.json').read_text())
matches=[j for j in jobs if j['job']==a.job]
if len(matches)!=1:raise SystemExit('Job ID is absent or ambiguous')
j=matches[0]
print(json.dumps(j,indent=2))
if not a.execute:raise SystemExit(0)
if sys.platform!='linux':raise SystemExit('Execute only in the prepared Linux cgroup scope')
cg=next(x.split(':',2)[2] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
cgpath=Path('/sys/fs/cgroup')/cg.lstrip('/')
if not cgpath.name.startswith('terreus-gates-'):raise SystemExit('Not inside the dedicated terreus-gates scope')
mem=(cgpath/'memory.max').read_text().strip();swap=(cgpath/'memory.swap.max').read_text().strip()
quota,period=(cgpath/'cpu.max').read_text().split()
if mem=='max' or int(mem)>8*1024**3 or swap!='0' or quota=='max' or int(quota)>int(period):raise SystemExit('Effective resource limits are not acceptable')
def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
input_path=j['argv'][j['argv'].index('--in')+1]
if sha(input_path)!=j['expected_file_sha256']:raise SystemExit('Input file changed')
lock=json.loads((OUT/'data/busco_lineage_lock.json').read_text())
lineage=Path(j['argv'][j['argv'].index('--lineage_dataset')+1])
for r in lock['files']:
 if sha(lineage/r['path'])!=r['sha256']:raise SystemExit('BUSCO lineage changed')
env=os.environ.copy();env.update(j['environment'])
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
version=subprocess.run(j['argv'][:2]+['--version'],env=env,capture_output=True,text=True,check=True).stdout.strip()
if version!='BUSCO 6.1.0':raise SystemExit('BUSCO version differs from prepared lock: '+version)
outparent=Path(j['argv'][j['argv'].index('--out_path')+1]);name=j['argv'][j['argv'].index('--out')+1]
if (outparent/name).exists():raise SystemExit('Output already exists; no overwrite or automatic restart')
outparent.mkdir(parents=True,exist_ok=True)
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
record={'job':j,'started_utc':stamp,'effective_cgroup':str(cgpath),'memory_max':mem,'cpu_max':[quota,period],'swap_max':swap,'version':version}
with (OUT/'logs'/f'{a.job}_{stamp}.log').open('w') as log:
 result=subprocess.run(j['argv'],env=env,stdout=log,stderr=subprocess.STDOUT)
record.update({'exit_code':result.returncode,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scientific_acceptance':'pending summary/output validation'})
(OUT/'logs'/f'{a.job}_{stamp}.json').write_text(json.dumps(record,indent=2))
raise SystemExit(result.returncode)
