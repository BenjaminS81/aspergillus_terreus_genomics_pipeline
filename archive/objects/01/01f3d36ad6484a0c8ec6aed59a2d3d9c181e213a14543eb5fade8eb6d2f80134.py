from pathlib import Path
import json,csv,collections
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11/lovastatin'
panel=json.loads((P/'reference_genes.json').read_text());core=[d['gene'] for d in panel if d['role']=='pathway_panel'];manifest=json.loads((P/'genome_manifest.json').read_text())
def unionlen(intervals):
    merged=[]
    for a,b in sorted(intervals):
        if merged and a<=merged[-1][1]+1:merged[-1][1]=max(b,merged[-1][1])
        else:merged.append([a,b])
    return sum(b-a+1 for a,b in merged)
def candidates(path,protein=False):
    hits=collections.defaultdict(list)
    for line in path.read_text().splitlines():
        a=line.split('\t');q,s=a[:2];pid=float(a[2]);length,qs,qe,ss,se=map(int,a[3:8]);ev,bs=map(float,a[8:10]);qlen,slen=map(int,a[10:12]);strand='+' if ss<=se else '-'
        hits[(q,s,strand)].append({'start':min(ss,se),'end':max(ss,se),'qs':qs,'qe':qe,'identity':pid,'length':length,'bitscore':bs,'qlen':qlen,'contig_length':slen})
    out=[]
    for (q,s,strand),hs in hits.items():
        groups=[]
        for h in sorted(hs,key=lambda h:h['start']):
            possible=[]
            for i,g in enumerate(groups):
                near=h['start']<=max(x['end'] for x in g)+2500
                ov=unionlen([(max(h['qs'],x['qs']),min(h['qe'],x['qe'])) for x in g if min(h['qe'],x['qe'])>=max(h['qs'],x['qs'])])
                if near and ov<.25*(h['qe']-h['qs']+1):possible.append(i)
            if possible:groups[possible[-1]].append(h)
            else:groups.append([h])
        for g in groups:
            cov=unionlen([(h['qs'],h['qe']) for h in g])/g[0]['qlen'];pid=sum(h['identity']*h['length'] for h in g)/sum(h['length'] for h in g)
            out.append({'gene':q,'contig':s,'strand':strand,'start':min(h['start'] for h in g),'end':max(h['end'] for h in g),'query_coverage':round(cov,4),'identity_percent':round(pid,3),'bitscore':round(sum(h['bitscore'] for h in g),2),'contig_length':g[0]['contig_length'],'hsp_count':len(g),'method':'tblastn' if protein else 'blastn','strong':cov>=.7 and pid>=(45 if protein else 75),'near_complete':cov>=.9 and pid>=(50 if protein else 80)})
    return out
rows=[];allcandidates={};rescue=[]
for m in manifest:
    sample=m['representative'];cs=candidates(P/'blastn'/(sample+'.tsv'))
    rescuepath=P/'tblastn'/(sample+'.tsv')
    if rescuepath.exists():
        aa=candidates(rescuepath,True)
        # Keep raw evidence separate; genomic candidates take precedence where strong.
        for c in aa:
            overlaps=[x for x in cs if x['gene']==c['gene'] and x['contig']==c['contig'] and x['strand']==c['strand'] and min(x['end'],c['end'])>=max(x['start'],c['start'])]
            if any(x['strong'] for x in overlaps):continue
            cs.append(c)
    allcandidates[sample]=cs
    best={g:max((c for c in cs if c['gene']==g),key=lambda c:(c['strong'],c['bitscore']),default=None) for g in core}
    missing=[g for g in core if not best[g] or not best[g]['strong']]
    if missing:rescue.append({**m,'genes':missing})
    # Every anchor/window is evaluated; clusters require colocalized evidence, not a list of unrelated best hits.
    anchors=[c for c in cs if c['gene'] in ['lovB','lovF'] and c['strong']];blocks=[]
    for anchor in anchors:
        window=[c for c in cs if c['gene'] in core and c['contig']==anchor['contig'] and abs((c['start']+c['end']-anchor['start']-anchor['end'])/2)<65000]
        chosen={g:max((c for c in window if c['gene']==g),key=lambda c:(c['strong'],c['bitscore']),default=None) for g in core}
        count=sum(bool(c and c['strong']) for c in chosen.values());score=sum(c['bitscore'] for c in chosen.values() if c and c['strong']);blocks.append((count,score,chosen,anchor))
    block=max(blocks,key=lambda x:(x[0],x[1])) if blocks else None
    if block and block[0]>=3:
        chosen=block[2];valid=[c for c in chosen.values() if c and c['strong']];low=min(c['start'] for c in valid);high=max(c['end'] for c in valid);contig=block[3]['contig'];length=block[3]['contig_length'];refB=chosen.get('lovB');orientation=1 if not refB or refB['strand']=='+' else -1
        if not refB or not refB['strong']:orientation=1 if chosen['lovF']['strand']=='+' else -1
        flank=[c for c in cs if c['gene'] not in core and c['strong'] and c['contig']==contig and c['end']>=low-20000 and c['start']<=high+20000]
        # Multiple alternative hits remain in all_candidates.json; diagrams use one best hit per reference flank.
        flank=list({g:max([c for c in flank if c['gene']==g],key=lambda c:c['bitscore']) for g in {c['gene'] for c in flank}}.values())
        row={'representative':sample,'samples':m['samples'],'subject':m['subject'],'status':'nine_gene_block' if block[0]==9 else 'partial_or_divergent_block','strong_colocalized_genes':block[0],'contig':contig,'start':low,'end':high,'contig_length':length,'orientation':orientation,'left_flank_available':min(20000,low-1),'right_flank_available':min(20000,length-high),'genes':chosen,'flanking_genes':flank}
    else:row={'representative':sample,'samples':m['samples'],'subject':m['subject'],'status':'no_convincing_colocalized_block','strong_colocalized_genes':0,'contig':None,'genes':{},'flanking_genes':[]}
    row['genome_wide_best_hits']=best;row['genome_wide_strong_genes']=sum(bool(c and c['strong']) for c in best.values());rows.append(row)
(P/'all_candidates.json').write_text(json.dumps(allcandidates,indent=2));(P/'cluster_screen.json').write_text(json.dumps(rows,indent=2));(P/'protein_rescue_manifest.json').write_text(json.dumps(rescue,indent=2))
print('BLOCKS',collections.Counter(x['status'] for x in rows),'RESCUE',len(rescue),'unique genomes',flush=True)
for x in rows:
    if x['status']!='nine_gene_block':print(x['representative'],x['status'],x['strong_colocalized_genes'],x['genome_wide_strong_genes'])
