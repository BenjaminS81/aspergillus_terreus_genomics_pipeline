from pathlib import Path
import json,hashlib,re,collections
BASE=Path(__file__).resolve().parent; P=BASE/'comparative_visualization/2026-11-11'; OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis'; ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
for d in ['sources','analysis','figures','lovastatin','provenance']:(P/d).mkdir(parents=True,exist_ok=True)
def fasta(p):
    name=None;s=[]
    for line in p.open():
        if line.startswith('>'):
            if name:yield name,''.join(s).upper()
            name=line[1:].split()[0];s=[]
        else:s.append(line.strip())
    if name:yield name,''.join(s).upper()
ref=ROOT/'data/reference/GCA_000149615.1_ASM14961v1_genomic.fna';gff=ROOT/'data/reference/genomic.gff'
genes=[]
for line in gff.open():
    a=line.rstrip().split('\t')
    if len(a)!=9 or a[2]!='gene':continue
    at=dict(x.split('=',1) for x in a[8].split(';') if '=' in x)
    if a[0]=='CH476609.1':genes.append({'id':at['locus_tag'],'contig':a[0],'start':int(a[3]),'end':int(a[4]),'strand':a[6]})
names=dict(zip(['ATEG_'+str(i).zfill(5) for i in range(9960,9969)],['lovA','lovB','lovG','lovC','lovD','lvrA','lovE','lovI','lovF']))
core=[x for x in genes if x['id'] in names];low=min(x['start'] for x in core);high=max(x['end'] for x in core)
selected=[x for x in genes if x['end']>=low-20000 and x['start']<=high+20000]
chrom=next(seq for n,seq in fasta(ref) if n=='CH476609.1')
nt=[]
for x in selected:
    x['gene']=names.get(x['id'],x['id']);x['role']='pathway_panel' if x['id'] in names else 'flanking_reference_gene'
    seq=chrom[x['start']-1:x['end']]
    if x['strand']=='-':seq=seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
    x['query_length']=len(seq);nt.append('>'+x['gene']+'\n'+seq+'\n')
(P/'lovastatin/reference_genomic_queries.fasta').write_text(''.join(nt))
(P/'lovastatin/reference_genes.json').write_text(json.dumps(selected,indent=2))
proteins=json.loads((P/'sources/uniprot_lov_cluster_nine.json').read_text())['results'];aa=[]
for x in proteins:
    locus=next(y['value'] for g in x['genes'] for y in g.get('orfNames',[]) if y['value'] in names)
    aa.append('>'+names[locus]+'\n'+x['sequence']['value']+'\n')
(P/'lovastatin/reference_protein_queries.fasta').write_text(''.join(aa))
sources=json.loads((OLD/'genome_sources.json').read_text());assets=json.loads((OLD/'assembly_variant_inventory.json').read_text());hashes={x['subject']:x.get('assembly_sequence_sha256') for x in assets};unique={}
for x in sources:
    h=hashes[x['subject']];assert h
    if h not in unique:unique[h]={'representative':x['sample'],'subject':x['subject'],'sequence_sha256':h,'samples':[]}
    unique[h]['samples'].append(x['sample'])
(P/'lovastatin/genome_manifest.json').write_text(json.dumps(list(unique.values()),indent=2))
(P/'sources/reference_locus_provenance.json').write_text(json.dumps({'reference_genome':str(ref),'genome_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'reference_gff':str(gff),'gff_sha256':hashlib.sha256(gff.read_bytes()).hexdigest(),'reference_contig':'CH476609.1','pathway_panel_start':low,'pathway_panel_end':high,'flank_bp':20000,'panel_definition':'Nine commonly compared pathway genes; this is not a claim that all regulatory/biosynthetic boundary genes are known. lovD/lvrA/lovI names mapped by the conserved ATEG locus series and the published nine-gene diagram, with UniProt proteins retained.','literature':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9814666/','mibig_limitation':'BGC0000088.4 has questionable quality and unknown completeness; lovB is not represented in that record.'},indent=2))
print('READY',len(unique),'distinct canonical genomes,',len(sources),'labels;',len(selected),'reference genes;',low,high)
