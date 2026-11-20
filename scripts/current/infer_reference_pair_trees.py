from pathlib import Path
import json,subprocess
P=Path(__file__).resolve().parent/'comparative_visualization/2026-11-11/reference_pair_combinations';runs=[]
for d in json.loads((P/'design.json').read_text()):
    name=d['name'];prefix=P/name;cmd=['/home/benjamin/miniforge3/envs/klebsiella_recombination/bin/iqtree2','-s',str(prefix)+'.fasta','-p',str(prefix)+'.nex','-pre',str(prefix),'-m','MFP+MERGE','-nt','1','-seed','20261111','-B','1000','-alrt','1000','-bnni','-o','REF_436_Aspergillus_flavipes_NRRL302']
    print('REFERENCE PAIR START',name,flush=True)
    if not (P/(name+'.iqtree')).exists():
        with (P/(name+'.stdout')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    assert 'Date and time:' in (P/(name+'.iqtree')).read_text();runs.append({'name':name,'command':cmd,'completed':True});(P/'execution.json').write_text(json.dumps(runs,indent=2));print('REFERENCE PAIR COMPLETE',name,flush=True)
