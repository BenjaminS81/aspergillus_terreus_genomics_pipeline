from pathlib import Path
import json,subprocess,time
P=Path(__file__).resolve().parent/'comparative_visualization/2026-11-11/lovastatin';out=P/'blastn';out.mkdir(exist_ok=True)
manifest=json.loads((P/'genome_manifest.json').read_text());runs=[]
for i,d in enumerate(manifest):
    target=out/(d['representative']+'.tsv');log=out/(d['representative']+'.log');done=out/(d['representative']+'.done')
    cmd=['/usr/bin/blastn','-task','blastn','-query',str(P/'reference_genomic_queries.fasta'),'-subject',d['subject'],'-evalue','1e-10','-dust','no','-perc_identity','60','-num_threads','1','-max_target_seqs','10000','-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen','-out',str(target)]
    if not done.exists():
        with log.open('w') as f:subprocess.run(cmd,stdout=f,stderr=f,check=True)
        done.write_text('completed\n')
    runs.append({'representative':d['representative'],'command':cmd,'completed':True})
    if i%6==0:print('BLAST',i+1,'/',len(manifest),d['representative'],flush=True)
(P/'blastn_execution.json').write_text(json.dumps(runs,indent=2))
print('ALL CANONICAL GENOME SCREENS COMPLETE',flush=True)
