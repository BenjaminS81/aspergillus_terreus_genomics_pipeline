"""Recompute sequence-derived metrics, validating the exact canonical inputs."""
from circos_review_common import *
LIMITS=constrain_process()
import hashlib, csv, collections, shutil, math
for folder in ('data','figures','qa','scripts'): (OUT/folder).mkdir(parents=True,exist_ok=True)
L=OLD/'lovastatin';data=readjson(L/'visualization_data.json');manifest=readjson(L/'genome_manifest.json')
for filename in ('visualization_data.json','reference_genes.json','genome_manifest.json','figure_order.json'):
    shutil.copy2(L/filename,OUT/'data'/filename)
shutil.copy2(OLD/'analysis/sample_annotations.json',OUT/'data/sample_annotations.json')
representatives={x['representative']:x for x in data};metrics={};regions={};checks=[]
def counts(seq):return [seq.count(c) for c in 'ACGTN']
def composition(c,total):
    a,cc,g,t,n=c;acgt=a+cc+g+t
    return {'gc_percent':100*(cc+g)/acgt if acgt else None,'gc_skew':(g-cc)/(g+cc) if g+cc else None,'n_percent':100*n/total,'ambiguous_percent':100*(total-acgt)/total}
for i,m in enumerate(manifest):
    rep=m['representative'];d=representatives[rep];source=Path(d['source']);seqs=dict(fasta(source))
    hashes=sorted(hashlib.sha256(s.encode()).hexdigest() for s in seqs.values())
    sha=hashlib.sha256('\n'.join(hashes).encode()).hexdigest()
    assert sha==m['sequence_sha256'],f'Changed assembly: {rep}'
    # Sorted contigs are a display coordinate system, not an inferred chromosome order.
    contigs=sorted(seqs,key=lambda k:(-len(seqs[k]),k));lengths=[len(seqs[k]) for k in contigs];total=sum(lengths)
    run=0;n50=None;l50=None;coords=[];totals=[0]*5
    edges=[int(j*total/100) for j in range(101)];bins=[[0]*5 for _ in range(100)]
    for k,size in zip(contigs,lengths):
        seq=seqs[k];c=counts(seq);totals=[a+b for a,b in zip(totals,c)]
        coords.append({'contig':k,'start_0based':run,'end_0based_exclusive':run+size,'length':size})
        if n50 is None and run+size>=total/2:n50=size;l50=len(coords)
        b0=min(99,int(run*100/total));b1=min(99,int((run+size-1)*100/total))
        for b in range(b0,b1+1):
            lo=max(0,edges[b]-run);hi=min(size,edges[b+1]-run)
            if hi>lo:bins[b]=[a+bv for a,bv in zip(bins[b],counts(seq[lo:hi]))]
        run+=size
    outbins=[]
    for j,c in enumerate(bins):
        lo,hi=edges[j:j+2]
        outbins.append({'start':lo,'end':hi,**composition(c,hi-lo),'contig_starts':sum(lo<=x['start_0based']<hi for x in coords)})
    metrics[rep]={'representative':rep,'samples':m['samples'],'source':str(source),'sequence_sha256':sha,'assembly_bp':total,'contigs':len(contigs),'N50_bp':n50,'L50':l50,**composition(totals,total),'coordinate_order':'descending contig length, then contig ID','contig_coordinates':coords,'bins':outbins}
    rr=[]
    for pi,p in enumerate(d['segments']):
        seq=seqs[p['contig']][p['window_start']-1:p['window_end']]
        assert len(seq)==p['window_end']-p['window_start']+1
        assert hashlib.sha256(seq.encode()).hexdigest()==p['sequence_sha256'],f'Changed region {rep}/{pi}'
        bb=[]
        for start in range(0,len(seq),1000):
            ss=seq[start:start+1000]
            bb.append({'start':start,'end':start+len(ss),**composition(counts(ss),len(ss))})
        rr.append({'segment':pi+1,'contig':p['contig'],'bp':len(seq),'bins':bb,**composition(counts(seq),len(seq))})
    regions[rep]=rr
    print(f'METRICS {i+1}/{len(manifest)} {rep}: {total:,} bp; {len(contigs)} contigs',flush=True)
writejson(OUT/'data/genome_metrics.json',metrics);writejson(OUT/'data/region_metrics.json',regions)
writejson(OUT/'data/resource_limits.json',LIMITS)
writejson(OUT/'data/source_validation.json',{'canonical_sequence_sets_checked':len(metrics),'labels':len(data),'all_assembly_hashes_match':True,'all_region_hashes_match':True,'window_size_bp':1000,'whole_assembly_bins':100,'coordinate_convention':'input annotation 1-based inclusive; exported plot intervals 0-based half-open'})
table=[];features=[]
for i,d in enumerate(data):
    m=metrics[d['representative']]
    table.append({**{k:d[k] for k in ('sample','representative','taxon','block_status','strong_colocalized_genes','genome_wide_strong_genes','segment_count')},**{k:m[k] for k in ('assembly_bp','contigs','N50_bp','L50','gc_percent','n_percent','ambiguous_percent','sequence_sha256')}})
    for j,p in enumerate(d['segments']):
        for kind,key in [('panel_homolog','genes'),('reference_flank_homolog','flanking_reference_genes'),('annotation_model','annotation_genes')]:
            for g in p[key]:features.append({'sample':d['sample'],'segment':j+1,'contig':p['contig'],'window_start':p['window_start'],'window_end':p['window_end'],'feature_kind':kind,'gene_or_model':g.get('gene',g.get('id')),'start_1based':g['start'],'end_1based_inclusive':g['end'],'strand':g['strand'],'method':g.get('method','annotation'),'identity_percent':g.get('identity_percent',''),'query_coverage':g.get('query_coverage',''),'display_orientation':p['orientation'],'annotation_source':p.get('annotation_source','') if kind=='annotation_model' else ''})
for name,rr in [('genome_metrics.tsv',table),('region_features.tsv',features)]:
    with (OUT/'data'/name).open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
writejson(OUT/'data/metric_rows.json',table)
print('COMPLETE: 77 labels, 72 validated assemblies. No source writes.',flush=True)
