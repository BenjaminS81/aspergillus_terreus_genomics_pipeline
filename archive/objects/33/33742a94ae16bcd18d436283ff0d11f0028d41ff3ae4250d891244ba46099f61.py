from pathlib import Path
import json,subprocess,concurrent.futures,hashlib,time
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
O=P/'genome_blast';O.mkdir(exist_ok=True)
mapping=json.loads((P/'genome_sources.json').read_text())
def job(d):
    output=O/(d['sample']+'.tsv');log=O/(d['sample']+'.log')
    cmd=['/usr/bin/blastn','-query',str(P/'genome_search_queries.fasta'),'-subject',d['subject'],'-task','blastn','-evalue','1e-15','-dust','no','-perc_identity','70','-max_target_seqs','10000','-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen','-out',str(output)]
    stamp=O/(d['sample']+'.source.txt')
    if not output.exists() or not log.exists() or 'error' in log.read_text().lower() or not stamp.exists() or stamp.read_text()!=d['subject']:
        r=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE);log.write_text(r.stdout+r.stderr)
        if r.returncode:raise RuntimeError((d['sample'],r.stderr))
        stamp.write_text(d['subject'])
    d['command']=cmd;d['subject_sha256']=hashlib.file_digest(open(d['subject'],'rb'),'sha256').hexdigest();d['output']=str(output.relative_to(P));d['hits']=len(output.read_text().splitlines())
    return d
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
    for d in ex.map(job,mapping):results.append(d);print('GENOME',len(results),'/',len(mapping),d['sample'],'hits',d['hits'],flush=True)
(P/'genome_search_manifest.json').write_text(json.dumps(results,indent=2))
