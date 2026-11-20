from pathlib import Path
import json,subprocess
P=Path(__file__).resolve().parent/'comparative_visualization/2026-11-11/lovastatin';O=P/'tblastn';O.mkdir(exist_ok=True)
seq={};k=None
for line in (P/'reference_protein_queries.fasta').read_text().splitlines():
    if line.startswith('>'):k=line[1:];seq[k]=''
    else:seq[k]+=line.strip()
manifest=json.loads((P/'protein_rescue_manifest.json').read_text());runs=[]
for i,d in enumerate(manifest):
    s=d['representative'];q=O/(s+'.faa');q.write_text(''.join('>'+g+'\n'+seq[g]+'\n' for g in d['genes']));target=O/(s+'.tsv');done=O/(s+'.done')
    cmd=['/usr/bin/tblastn','-query',str(q),'-subject',d['subject'],'-evalue','1e-10','-seg','yes','-num_threads','1','-max_target_seqs','10000','-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen','-out',str(target)]
    if not done.exists():
        with (O/(s+'.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=f,check=True)
        done.write_text('completed\n')
    runs.append({'sample':s,'command':cmd,'completed':True});(O/'execution.json').write_text(json.dumps(runs,indent=2));print('PROTEIN RESCUE',i+1,'/',len(manifest),s,flush=True)
