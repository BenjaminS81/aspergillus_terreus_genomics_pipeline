from pathlib import Path
import json,collections,csv
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
def fasta(p):
    d={};n=None
    for ln in p.read_text().splitlines():
        if ln.startswith('>'):n=ln[1:].split()[0];d[n]=''
        else:d[n]+=ln.strip()
    return d
comp=str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN')
report=[]
for l in ['BenA','CaM','RPB2']:
    seqs=fasta(P/(l+'_local_ungapped.fasta'));groups=collections.defaultdict(list)
    for line in (P/(l+'_reference_to_local.blast.tsv')).read_text().splitlines():
        v=line.split('\t');d=dict(zip(['query','sample','identity','length','mismatch','gapopen','qs','qe','ss','se','evalue','bitscore','qlen','slen'],v))
        for k in ['length','qs','qe','ss','se','qlen','slen']:d[k]=int(d[k])
        d['bitscore']=float(d['bitscore']);d['identity']=float(d['identity']);d['strand']=1 if d['ss']<d['se'] else -1
        groups[(d['sample'],d['query'],d['strand'])].append(d)
    best={}
    for (sample,query,strand),hsps in groups.items():
        kept=[];covered=set()
        for h in sorted(hsps,key=lambda x:-x['bitscore']):
            region=set(range(min(h['ss'],h['se']),max(h['ss'],h['se'])+1))
            if len(region-covered)/len(region)<.8:continue
            kept.append(h);covered.update(region)
        score=sum(h['bitscore'] for h in kept)
        if sample not in best or score>best[sample][0]:best[sample]=(score,query,strand,kept)
    cropped={}
    for sample,seq in seqs.items():
        if sample not in best:report.append({'locus':l,'sample':sample,'status':'no_blast_hit'});continue
        score,query,strand,hsps=best[sample];lo=min(min(h['ss'],h['se']) for h in hsps);hi=max(max(h['ss'],h['se']) for h in hsps)
        qcovered=set(p for h in hsps for p in range(h['qs'],h['qe']+1));ident=sum(h['identity']*h['length'] for h in hsps)/sum(h['length'] for h in hsps)
        status='pass' if len(qcovered)>=120 and ident>=70 else 'hold_low_overlap'
        cut=seq[lo-1:hi]
        if strand==-1:cut=cut.translate(comp)[::-1]
        if status=='pass':cropped[sample]=cut
        report.append({'locus':l,'sample':sample,'status':status,'best_reference':query,'strand':strand,'start_1based':lo,'end_1based_inclusive':hi,'original_length':len(seq),'cropped_length':len(cut),'matched_reference_bases':len(qcovered),'reference_length':hsps[0]['qlen'],'weighted_identity_percent':round(ident,3),'hsp_count':len(hsps),'bitscore_sum':score})
    (P/(l+'_local_window.fasta')).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in cropped.items()))
    refs=(P/(l+'_references.fasta')).read_text();(P/(l+'_combined_unaligned.fasta')).write_text(refs+''.join('>'+k+'\n'+v+'\n' for k,v in cropped.items()))
    print(l,'passed',len(cropped),'ranges',min(map(len,cropped.values())),max(map(len,cropped.values())))
(P/'local_window_validation.json').write_text(json.dumps(report,indent=2))
cols=list(dict.fromkeys(k for d in report for k in d))
with (P/'local_window_validation.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(report)
