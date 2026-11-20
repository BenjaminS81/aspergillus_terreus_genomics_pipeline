"""Dependency launcher: no scientific job starts until the pilot service has exited."""
from pathlib import Path
import subprocess,time,os,sys,json
OUT=Path(__file__).resolve().parents[1]
PILOT='terreus-gates-pilot2-20261118.service'
print('Waiting for the single pilot; no analysis running in this dispatcher',flush=True)
while True:
 r=subprocess.run(['systemctl','--user','show',PILOT,'-p','ActiveState','--value'],capture_output=True,text=True)
 if r.stdout.strip() not in {'active','activating','deactivating','reloading'}:break
 time.sleep(15)
if not (OUT/'reports/BUSCO_genome_SRR9691936.validated.json').exists():
 (OUT/'reports/pipeline_launch_blocked.json').write_text(json.dumps({'reason':'Pilot did not produce a validated result; full queue not started','pilot_service':PILOT},indent=2))
 raise SystemExit(1)
print('Pilot passed and exited; starting sequential gate work',flush=True)
os.execv(sys.executable,[sys.executable,'-u',str(OUT/'scripts/worker.py'),'full'])
