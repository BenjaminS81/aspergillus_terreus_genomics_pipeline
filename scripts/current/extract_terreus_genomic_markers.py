from pathlib import Path
import json,collections,hashlib,csv,sys
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
OUT=P/(sys.argv[2] if len(sys.argv)>2 else 'genomic_markers');OUT.mkdir(exist_ok=True)
sources=json.loads((P/(sys.argv[1] if len(sys.argv)>1 else 'genome_search_manifest.json')).read_text());allhits=[];calls=[];fastas={l:{} for l in ['BenA','CaM','RPB2']}
comp=str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN')
for src in sources:
    groups=collections.defaultdict(list)
    for line in (P/src['output']).read_text().splitlines():
        v=line.split('\t');h=dict(zip(['query','contig','identity','length','mismatch','gapopen','qs','qe','ss','se','evalue','bitscore','qlen','slen'],v))
        for k in ['length','qs','qe','ss','se','qlen','slen']:h[k]=int(h[k])
        for k in ['identity','bitscore']:h[k]=float(h[k])
        h['locus']=h['query'].split('_')[0];h['strand']=1 if h['ss']<h['se'] else -1;h['sample']=src['sample'];h['start']=min(h['ss'],h['se']);h['end']=max(h['ss'],h['se'])
        groups[(h['query'],h['contig'],h['strand'])].append(h);allhits.append(h)
    candidates=collections.defaultdict(list)
    for (query,contig,strand),hs in groups.items():
        locus=query.split('_')[0]
        clusters=[]
        for h in sorted(hs,key=lambda x:x['start']):
            repeated_query=False
            if clusters:
                prior=set(p for x in clusters[-1] for p in range(x['qs'],x['qe']+1));q=set(range(h['qs'],h['qe']+1))
                repeated_query=len(prior&q)/len(q)>.5 and h['start']>max(x['end'] for x in clusters[-1])
            if not clusters or h['start']-max(x['end'] for x in clusters[-1])>2000 or repeated_query:clusters.append([])
            clusters[-1].append(h)
        for cluster in clusters:
            keep=[];covered=set()
            for h in sorted(cluster,key=lambda x:-x['bitscore']):
                q=set(range(h['qs'],h['qe']+1))
                if len(q-covered)/len(q)<.8:continue
                covered|=q;keep.append(h)
            n=sum(h['length'] for h in keep);identity=sum(h['identity']*h['length'] for h in keep)/n
            c={'sample':src['sample'],'locus':locus,'query':query,'contig':contig,'strand':strand,'start_1based':min(h['start'] for h in keep),'end_1based_inclusive':max(h['end'] for h in keep),'query_coverage_percent':round(100*len(covered)/keep[0]['qlen'],3),'weighted_identity_percent':round(identity,3),'bitscore_sum':sum(h['bitscore'] for h in keep),'hsp_count':len(keep)}
            c['strong_candidate']=c['query_coverage_percent']>=70 and identity>=85;c['near_complete']=c['query_coverage_percent']>=90 and identity>=85
            candidates[locus].append(c)
    wanted=set(c['contig'] for cs in candidates.values() for c in cs);seqs={};name=None
    with open(src['subject']) as f:
        for line in f:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in wanted:seqs[name]=[]
            elif name in wanted:seqs[name].append(line.strip().upper())
    seqs={k:''.join(v) for k,v in seqs.items()}
    for locus in ['BenA','CaM','RPB2']:
        cs=[]
        for c in sorted(candidates[locus],key=lambda c:(not c['strong_candidate'],-c['bitscore_sum'])):
            if any(c['contig']==x['contig'] and c['strand']==x['strand'] and min(c['end_1based_inclusive'],x['end_1based_inclusive'])>=max(c['start_1based'],x['start_1based']) for x in cs):continue
            cs.append(c)
        strong=[c for c in cs if c['strong_candidate']]
        status='single_strong_candidate' if len(strong)==1 and strong[0]['near_complete'] else 'partial_single_candidate' if len(strong)==1 else 'multiple_strong_candidates_hold' if len(strong)>1 else 'no_strong_candidate_hold'
        call={'sample':src['sample'],'locus':locus,'status':status,'strong_candidate_count':len(strong),'total_candidate_regions':len(cs),'in_core_tree':src['in_core_tree'],'genome_source':src['subject'],'source_sha256':src['subject_sha256'],'candidates':cs}
        if len(strong)==1:
            c=strong[0];seq=seqs[c['contig']][c['start_1based']-1:c['end_1based_inclusive']]
            if c['strand']==-1:seq=seq.translate(comp)[::-1]
            call.update({k:v for k,v in c.items() if k not in ['sample','locus']});call['sequence']=seq
            if set(seq)-set('ACGTRYKMSWBDHVN'):raise ValueError((src['sample'],locus,set(seq)))
            fastas[locus]['LOCAL_'+src['sample']]=seq
        calls.append(call)
    print('EXTRACT',src['sample'],[(c['locus'],c['status']) for c in calls[-3:]],flush=True)
for l,seqs in fastas.items():(OUT/(l+'_genomic.fasta')).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in seqs.items()))
(OUT/'locus_calls.json').write_text(json.dumps(calls,indent=2));(OUT/'all_blast_hits.json').write_text(json.dumps(allhits,indent=2))
summary=[{k:v for k,v in c.items() if k not in ['candidates','sequence']} for c in calls];cols=list(dict.fromkeys(k for d in summary for k in d))
with (OUT/'locus_calls.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(summary)
print('CALL SUMMARY',collections.Counter(c['status'] for c in calls))
