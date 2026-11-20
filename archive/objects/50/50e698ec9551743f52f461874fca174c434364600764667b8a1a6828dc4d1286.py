"""Verify a user-pause checkpoint; resume ONLY with explicit --resume."""
from pathlib import Path
import argparse,datetime,hashlib,json,subprocess,shutil
OUT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group();g.add_argument('--check',action='store_true');g.add_argument('--resume',action='store_true');args=parser.parse_args()
 flag=OUT/'PAUSED_BY_USER.json';pause=load(flag);folder=Path(pause['checkpoint']);snap=load(folder/'checkpoint.json');unit=snap['unit']
 assert unit=='terreus-gates-all20gb-20261118.service'
 for r in snap['validated_results']:
  p=Path(r['validation_file']);assert sha(p)==r['sha256'],str(p)
  result=load(p);summary=Path(result['summary']);assert sha(summary)==result['summary_sha256']
  full_table=summary.parent/'run_eurotiales_odb12.2/full_table.tsv';assert sha(full_table)==result['full_table_sha256']
  target=folder/'validated_results'/p.name;target.parent.mkdir(exist_ok=True)
  if target.exists():assert sha(target)==r['sha256']
  else:shutil.copy2(p,target)
 for r in snap['busco_disk_checkpoints']:
  assert sha(Path(r['saved']))==r['sha256'];assert sha(Path(r['original']))==r['sha256']
 props=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',unit,'-p','FreezerState','-p','MainPID','-p','ControlGroup','-p','MemoryMax','-p','MemorySwapMax','-p','CPUQuotaPerSecUSec'],text=True).splitlines())
 live=props.get('FreezerState')=='frozen' and props.get('MainPID')==snap['service']['MainPID'] and props.get('ControlGroup')==snap['service']['ControlGroup']
 checked={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'validated_completed_genomes':len(snap['validated_results']),'summaries_and_full_tables_hash_verified':True,'saved_busco_checkpoints':len(snap['busco_disk_checkpoints']),'live_resume_available':live,'service':props,'resumed':False}
 (folder/'verification.json').write_text(json.dumps(checked,indent=2))
 if not args.resume:print(json.dumps(checked,indent=2));return
 assert live,'Live state unavailable: follow cold recovery notes; do not start duplicate jobs'
 assert int(props['MemoryMax'])<=20_000_000_000 and props['MemorySwapMax']=='0' and props['CPUQuotaPerSecUSec']=='1s'
 available=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
 assert available>=20_000_000_000+12*1024**3,'Insufficient RAM headroom for safe resumption'
 current=load(OUT/'status.json');original=load(folder/'status.json')
 original['resumed_at']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 (OUT/'status.json').write_text(json.dumps(original,indent=2))
 try:subprocess.run(['systemctl','--user','thaw',unit],check=True)
 except Exception:
  (OUT/'status.json').write_text(json.dumps(current,indent=2));raise
 flag.rename(folder/'PAUSE_RESOLVED.json')
 checked.update({'resumed':True,'resumed_at':original['resumed_at']});(folder/'verification.json').write_text(json.dumps(checked,indent=2))
 print('Resumed the existing Terreus service; no additional worker launched.')
if __name__=='__main__':main()
