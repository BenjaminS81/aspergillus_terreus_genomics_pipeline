from pathlib import Path
import json,subprocess,time
P=Path(__file__).resolve().parent/'comparative_visualization/2026-11-11/core66_combinations';exe='/home/benjamin/miniforge3/envs/klebsiella_recombination/bin/iqtree2';runs=[]
for d in json.loads((P/'design.json').read_text())['combinations']:
    name=d['name'];prefix=P/name
    cmd=[exe,'-s',str(prefix)+'.fasta','-pre',str(prefix),'-nt','1','-seed','20261111','-B','1000','-alrt','1000','-bnni','-o','LOCAL_SRR10507650']
    if len(d['genes'])>1:cmd+=['-p',str(prefix)+'.nex','-m','MFP+MERGE']
    else:cmd+=['-m','MFP','-mem','1G']
    print('START',name,flush=True);t=time.time()
    if not (P/(name+'.iqtree')).exists():
        with (P/(name+'.stdout')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    assert 'Date and time:' in (P/(name+'.iqtree')).read_text()
    runs.append({'name':name,'command':cmd,'completed':True,'elapsed_seconds_this_invocation':round(time.time()-t,2)})
    (P/'execution.json').write_text(json.dumps(runs,indent=2));print('COMPLETE',name,flush=True)
