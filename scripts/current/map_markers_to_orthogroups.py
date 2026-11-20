from pathlib import Path
import json,csv,re,collections,hashlib
BASE=Path(__file__).resolve().parent;P=BASE/'taxonomy_reference_build/2026-11-10/analysis';ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
calls=json.loads((P/'genomic_markers/locus_calls.json').read_text());comp=str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN')
def fasta(path):
    name=None;seq=[]
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(seq).upper()
                name=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
    if name is not None:yield name,''.join(seq).upper()
def geneid(s):return re.sub(r'-T\d+$','',s.split('|')[-1])
old=[]
for locus,fn in [('BenA','BENA56.fa'),('CaM','CAM1.fa'),('RPB2','RPB2.fa')]:
    for name,seq in fasta(ROOT/'results/Houbraken_genes'/fn):
        sample,gene=name.split('|',1);sample=re.sub(r'^(terreus|hortae|pseudoterreus|alabamensis|aureoterreus|floccosus)_','',sample)
        if gene!='MISSING':old.append({'sample':sample,'locus':locus,'old_gene':gene,'original_header':name,'sequence_length':len(seq)})
mapped=[];targets=collections.defaultdict(set)
for d in old:targets[d['sample']].add(geneid(d['old_gene']))
for sample in sorted({d['sample'] for d in calls}):
    c=[d for d in calls if d['sample']==sample and d.get('sequence')]
    folders=list((ROOT/'results/0_2_repeat_mask_funannotate').glob('*/'+sample+'/*/annotate_results'))
    if not folders:
        for d in c:mapped.append({'sample':sample,'locus':d['locus'],'status':'annotation_not_found'})
        continue
    folder=sorted(folders)[0];gffs=list(folder.glob('*.gff3'));genomes=list(folder.glob('*.fasta'))
    if len(gffs)!=1 or len(genomes)!=1:raise ValueError((sample,gffs,genomes))
    genes=[];mrna={}
    for line in gffs[0].read_text().splitlines():
        if not line or line.startswith('#'):continue
        a=line.split('\t')
        if len(a)!=9 or a[2] not in ['gene','mRNA']:continue
        attr=dict(item.split('=',1) for item in a[8].split(';') if '=' in item)
        d={'contig':a[0],'start':int(a[3]),'end':int(a[4]),'strand':a[6],'attributes':attr}
        if a[2]=='gene':genes.append(d)
        else:mrna[attr.get('Parent','')]=attr
    hits=collections.defaultdict(list)
    for name,seq in fasta(genomes[0]):
        for d in c:
            for strand,query in [('+',d['sequence']),('-',d['sequence'].translate(comp)[::-1])]:
                start=seq.find(query)
                while start>=0:
                    hits[d['locus']].append((name,start+1,start+len(query),strand));start=seq.find(query,start+1)
    for d in c:
        hh=hits[d['locus']];matches=[]
        for name,start,end,strand in hh:
            for g in genes:
                ov=min(end,g['end'])-max(start,g['start'])+1
                if g['contig']==name and g['strand']==strand and ov>0:
                    ident=g['attributes']['ID'];targets[sample].add(geneid(ident));matches.append({'gene_id':ident,'marker_contig':name,'marker_start':start,'marker_end':end,'gene_start':g['start'],'gene_end':g['end'],'overlap_bp':ov,'gene_name':g['attributes'].get('Name',''),'mrna_annotation':mrna.get(ident,{})})
        mapped.append({'sample':sample,'locus':d['locus'],'status':'mapped' if len(matches)==1 else 'multiple_gene_models' if matches else 'no_exact_annotation_mapping','annotation_gff':str(gffs[0]),'marker_exact_matches':len(hh),'genes':matches})
    # Preserve annotations of the old extracts for direct comparison.
    for d in old:
        if d['sample']==sample:
            g=next((g for g in genes if g['attributes']['ID']==geneid(d['old_gene'])),None)
            if g:d['annotation']={**g,'mrna':mrna.get(g['attributes']['ID'],{})}
    print('ANNOTATION MAP',sample,flush=True)
orth=ROOT/'results/0_3_orthofinder/terrei3/Results_Sep28/Orthogroups/Orthogroups.tsv';found={};ogrows={}
with orth.open(newline='') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        for sample,genes in targets.items():
            cell=row.get(sample+'_transcripts','');members=[x.strip() for x in cell.split(',') if x.strip()]
            matched={geneid(x) for x in members}&genes
            if matched:
                ogrows[row['Orthogroup']]=row
                for gid in matched:found[(sample,gid)]=row['Orthogroup']
for d in old:d['orthogroup']=found.get((d['sample'],geneid(d['old_gene'])))
for d in mapped:
    for g in d.get('genes',[]):g['orthogroup']=found.get((d['sample'],geneid(g['gene_id'])))
out=P/'orthology_crosscheck';out.mkdir(exist_ok=True)
(out/'old_marker_genes.json').write_text(json.dumps(old,indent=2));(out/'validated_marker_gene_mapping.json').write_text(json.dumps(mapped,indent=2));(out/'relevant_orthogroup_rows.json').write_text(json.dumps(ogrows,indent=2))
summary={l:{'old':dict(collections.Counter(d.get('orthogroup') for d in old if d['locus']==l)),'validated':dict(collections.Counter(g.get('orthogroup') for d in mapped if d['locus']==l for g in d.get('genes',[])))} for l in ['BenA','CaM','RPB2']}
(out/'orthogroup_summary.json').write_text(json.dumps(summary,indent=2));print('ORTHOGROUPS',json.dumps(summary))
print('NIH2624 OLD',json.dumps([d for d in old if d['sample']=='NIH2624']))
print('NIH2624 VALIDATED',json.dumps([d for d in mapped if d['sample']=='NIH2624']))
