from pathlib import Path
import json,hashlib,collections,subprocess
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
assets=json.loads((P/'assembly_asset_inventory.json').read_text());canonical=json.loads((P/'genome_search_manifest.json').read_text());cachepath=P/'assembly_sequence_hash_cache.json';cache=json.loads(cachepath.read_text()) if cachepath.exists() else {}
def digest(path):
    h=hashlib.sha256();seqhashes=[];count=0;length=0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if count:seqhashes.append(h.hexdigest());h=hashlib.sha256()
                count+=1
            else:
                s=line.strip().upper();h.update(s.encode());length+=len(s)
    if count:seqhashes.append(h.hexdigest())
    return {'assembly_sequence_sha256':hashlib.sha256('\n'.join(sorted(seqhashes)).encode()).hexdigest(),'contigs':count,'bases':length}
for i,d in enumerate(assets):
    if not d['bytes']:d['audit_status']='empty_file';continue
    p=Path(d['subject']);key=str(p);st=p.stat();tag=str(st.st_size)+':'+str(st.st_mtime_ns)
    if key not in cache or cache[key]['tag']!=tag:cache[key]={'tag':tag,**digest(p)}
    d.update(cache[key]);d['audit_status']='hashed'
    if i%25==0:cachepath.write_text(json.dumps(cache,indent=2));print('ASSEMBLY FILES',i+1,'/',len(assets),flush=True)
cachepath.write_text(json.dumps(cache,indent=2));(P/'assembly_variant_inventory.json').write_text(json.dumps(assets,indent=2))
bydigest=collections.defaultdict(list)
for d in assets:
    if d.get('assembly_sequence_sha256'):bydigest[d['assembly_sequence_sha256']].append(d)
known={cache[d['subject']]['assembly_sequence_sha256']:d for d in canonical}
jobs=[]
for sha,ds in bydigest.items():
    if sha in known:continue
    chosen=sorted(ds,key=lambda d:(d['priority'],d['subject']))[0]
    jobs.append({'sample':'VAR_'+sha[:12],'subject':chosen['subject'],'assembly_sequence_sha256':sha,'sample_aliases':sorted({d['sample'] for d in ds}),'in_core_tree':False})
V=P/'variant_blast';V.mkdir(exist_ok=True)
for i,d in enumerate(jobs):
    output=V/(d['sample']+'.tsv');log=V/(d['sample']+'.log')
    cmd=['/usr/bin/blastn','-query',str(P/'genome_search_queries.fasta'),'-subject',d['subject'],'-task','blastn','-evalue','1e-15','-dust','no','-perc_identity','70','-max_target_seqs','10000','-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen','-out',str(output)]
    if not output.exists() or not log.exists():
        r=subprocess.run(cmd,text=True,capture_output=True);log.write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    d.update({'output':str(output.relative_to(P)),'subject_sha256':hashlib.file_digest(open(d['subject'],'rb'),'sha256').hexdigest(),'command':cmd})
    print('ADDITIONAL ASSEMBLY VARIANT',i+1,'/',len(jobs),d['sample_aliases'],flush=True)
(P/'variant_search_manifest.json').write_text(json.dumps(jobs,indent=2));(P/'assembly_sequence_groups.json').write_text(json.dumps({sha:{'samples':sorted({d['sample'] for d in ds}),'files':[d['subject'] for d in ds],'searched_as':known[sha]['sample'] if sha in known else next(j['sample'] for j in jobs if j['assembly_sequence_sha256']==sha)} for sha,ds in bydigest.items()},indent=2))
print('UNIQUE ASSEMBLY SEQUENCE SETS',len(bydigest),'CANONICAL SEQUENCE SETS',len(known),'ADDITIONAL',len(jobs))
