from pathlib import Path
import json,collections,csv,re,hashlib
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';L=P/'lovastatin';OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis';ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
rows=json.loads((L/'cluster_screen.json').read_text());candidates=json.loads((L/'all_candidates.json').read_text());panel=json.loads((L/'reference_genes.json').read_text());ref={x['gene']:x for x in panel};core=[x['gene'] for x in panel if x['role']=='pathway_panel'];ann=json.loads((P/'analysis/sample_annotations.json').read_text());maps=json.loads((OLD/'orthology_crosscheck/validated_marker_gene_mapping.json').read_text());calls=json.loads((OLD/'genomic_markers/locus_calls.json').read_text());source={x['sample']:x for x in json.loads((OLD/'genome_sources.json').read_text())}
def winpath(s):return Path(s.replace('/mnt/d/','D:/'))
def fasta(p):
    k=None;s=[]
    for line in p.open():
        if line.startswith('>'):
            if k:yield k,''.join(s).upper()
            k=line[1:].split()[0];s=[]
        else:s.append(line.strip())
    if k:yield k,''.join(s).upper()
regions=[];summary=[];gene_rows=[];seqdir=L/'region_fastas';seqdir.mkdir(exist_ok=True)
for row in rows:
    rep=row['representative'];best=row['genome_wide_best_hits'];valid=[c for c in best.values() if c and c['strong']]
    # The selected primary block remains preferred over a stronger off-cluster paralog.
    for g,c in row['genes'].items():
        if c and c['strong']:valid=[x for x in valid if x['gene']!=g]+[c]
    bycontig=collections.defaultdict(list)
    for c in valid:bycontig[c['contig']].append(c)
    pieces=[]
    for contig,genes in bycontig.items():
        blocks=[]
        for g in sorted(genes,key=lambda x:x['start']):
            if blocks and g['start']-max(x['end'] for x in blocks[-1])<=25000:blocks[-1].append(g)
            else:blocks.append([g])
        for block in blocks:
            anchor=next((x for x in block if x['gene']=='lovB'),next((x for x in block if x['gene']=='lovF'),block[0]));orient=1 if anchor['strand']==ref[anchor['gene']]['strand'] else -1
            low=min(x['start'] for x in block);high=max(x['end'] for x in block);size=anchor['contig_length'];lo=max(1,low-20000);hi=min(size,high+20000)
            flank=[c for c in candidates[rep] if c['gene'] not in core and c['strong'] and c['contig']==contig and c['end']>=lo and c['start']<=hi];flank=list({g:max([c for c in flank if c['gene']==g],key=lambda x:x['bitscore']) for g in {c['gene'] for c in flank}}.values())
            pieces.append({'contig':contig,'start':low,'end':high,'window_start':lo,'window_end':hi,'contig_length':size,'orientation':orient,'anchor':anchor['gene'],'genes':block,'flanking_reference_genes':flank,'annotation_genes':[],'context_status':'reference-homology flanks only','fragment_edge':low<=2000 or size-high<=2000,'reference_offset':ref[anchor['gene']]['start']-(anchor['start'] if orient==1 else -anchor['end'])})
    pieces.sort(key=lambda x:(-len(x['genes']),x['contig'],x['start']))
    # Use native GFF beside the public assembly, or an annotation demonstrably in the same coordinate system.
    subject=winpath(row['subject']);gff=subject.parent/'genomic.gff'
    if not gff.exists():
        good=[]
        for m in maps:
            if m['sample']!=rep or 'annotation_gff' not in m:continue
            c=next((c for c in calls if c['sample']==rep and c['locus']==m['locus'] and c.get('sequence')),None)
            if c and any(g['marker_contig']==c['contig'] and g['marker_start']==c['start_1based'] and g['marker_end']==c['end_1based_inclusive'] for g in m.get('genes',[])):good.append(Path(m['annotation_gff']))
        gff=good[0] if good else None
    if gff and gff.exists():
        for line in gff.open():
            a=line.rstrip().split('\t')
            if len(a)!=9 or a[2]!='gene':continue
            at=dict(x.split('=',1) for x in a[8].split(';') if '=' in x);start,end=int(a[3]),int(a[4])
            for piece in pieces:
                if a[0]==piece['contig'] and end>=piece['window_start'] and start<=piece['window_end']:
                    piece['annotation_genes'].append({'id':at.get('locus_tag',at.get('ID','')),'start':start,'end':end,'strand':a[6]});piece['context_status']='matching-coordinate annotated genes plus reference homologs';piece['annotation_source']=str(gff)
    wanted={x['contig'] for x in pieces};seqs={k:s for k,s in fasta(subject) if k in wanted};out=[]
    if gff and gff.exists() and gff.name!='genomic.gff':
        annotation_fastas=list(gff.parent.glob('*.fasta'))
        annotation_seqs={k:s for k,s in fasta(annotation_fastas[0]) if k in wanted} if len(annotation_fastas)==1 else {}
        for piece in pieces:
            a,b=piece['window_start']-1,piece['window_end'];same=annotation_seqs.get(piece['contig'],'')[a:b]==seqs[piece['contig']][a:b]
            piece['annotation_window_exact_match']=same
            if not same:piece['annotation_genes']=[];piece['context_status']='reference-homology flanks only; annotation-coordinate check failed'
    for i,piece in enumerate(pieces):
        seq=seqs[piece['contig']][piece['window_start']-1:piece['window_end']];assert len(seq)==piece['window_end']-piece['window_start']+1
        out.append(f'>{rep}|segment{i+1}|{piece["contig"]}:{piece["window_start"]}-{piece["window_end"]}|genomic_orientation\n'+seq+'\n');piece['sequence_sha256']=hashlib.sha256(seq.encode()).hexdigest()
    if out:(seqdir/(rep+'.fasta')).write_text(''.join(out))
    for sample in row['samples']:
        d={'sample':sample,'representative':rep,'taxon':ann['LOCAL_'+sample]['taxon'],'source':str(subject),'block_status':row['status'],'strong_colocalized_genes':row['strong_colocalized_genes'],'genome_wide_strong_genes':row['genome_wide_strong_genes'],'segment_count':len(pieces),'segments':pieces,'genome_wide_best_hits':best};regions.append(d)
        summary.append({k:v for k,v in d.items() if k not in ['segments','genome_wide_best_hits']})
        for g in core:
            c=next((c for p in pieces for c in p['genes'] if c['gene']==g),None)
            gene_rows.append({'sample':sample,'gene':g,'status':'strong_homolog' if c else 'below_threshold_or_not_detected','cluster_status':row['status'],'contig':c['contig'] if c else '', 'start':c['start'] if c else '', 'end':c['end'] if c else '', 'strand':c['strand'] if c else '', 'query_coverage':c['query_coverage'] if c else '', 'identity_percent':c['identity_percent'] if c else '', 'method':c['method'] if c else ''})
    print('CONTEXT',rep,len(pieces),'segments',flush=True)
(L/'visualization_data.json').write_text(json.dumps(regions,indent=2));(L/'coverage_summary.json').write_text(json.dumps(summary,indent=2))
for fn,rr in [('coverage_summary.tsv',summary),('gene_evidence.tsv',gene_rows)]:
    with (L/fn).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
print('ALL LABELS',len(regions),collections.Counter(d['block_status'] for d in regions))
