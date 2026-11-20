"""Read-only scientific audit; writes derived evidence only to the C: workspace."""
from pathlib import Path
import json, csv, re, hashlib, collections, datetime, base64
import numpy as np

ROOT = Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
OUT = Path(__file__).resolve().parent / 'review_reports' / '2026-11-10'
OUT.mkdir(parents=True, exist_ok=True)
assert ROOT not in OUT.parents
E = {'reviewed_at':datetime.datetime.now().astimezone().isoformat(), 'source':str(ROOT), 'scope':'Read-only targeted content and code audit, not exhaustive validation of every intermediate or raw read.'}

def table(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def fasta(p):
    name=None; chunks=[]
    with Path(p).open(encoding='utf-8-sig') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(chunks)
                name=line[1:].strip();chunks=[]
            else:chunks.append(line.strip())
    if name is not None:yield name,''.join(chunks)

def tips(p):return re.findall(r'(?<=[(,])([^():,;]+)(?=:|[,);])',Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def norm(s):
    s=s.split('|')[0].removesuffix('_transcripts').removesuffix('.transcripts')
    for prefix in ['terreus_','hortae_','alabamensis_','aureoterreus_','floccosus_','pseudoterreus_']:
        if s.startswith(prefix):s=s[len(prefix):];break
    return s

# Exact current run evidence.
bundles={
 'genome':'results/0_4_mafft_iqtree/Sep28/concatenated.fa',
 'markers':'results/Houbraken_genes/merged/concatenated.fa',
 'rpb2':'results/Houbraken_genes/mafft_iqtree/RPB2_aligned.fa',
 'snp':'results/1_5_snp_iqtree/filtered.min0.phy'}
E['trees']={}
for key,rel in bundles.items():
    p=ROOT/rel; log=Path(str(p)+'.log').read_text(errors='replace')
    tree=Path(str(p)+'.treefile');names=tips(tree)
    E['trees'][key]={'path':str(tree),'tips':names,'normalized_tips':list(map(norm,names)),
      'sha256':sha(tree),'log_evidence':[l for l in log.splitlines() if any(s in l for s in ['IQ-TREE version','Command:','Seed:','Alignment has','parsimony-informative','TOTAL','Best-fit model','Total wall-clock','WARNING:'])]}
    mod=Path(str(p)+'.mod.nwk')
    if mod.exists():
        # Verify exact relabel-only transformation by removing labels from full Newick text.
        a=re.sub(r'(?<=[(,])([^():,;]+)(?=:|[,);])','TIP',tree.read_text()).strip()
        b=re.sub(r'(?<=[(,])([^():,;]+)(?=:|[,);])','TIP',mod.read_text()).strip()
        E['trees'][key]['modified_variant']={'path':str(mod),'tips':tips(mod),'identical_structure_branch_text':a==b,'same_normalized_tip_order':list(map(norm,names))==list(map(norm,tips(mod)))}

# Gene-count inventory; neither zero nor multiple-copy values are interpreted as biology here.
op=ROOT/'results/0_3_orthofinder/terrei3/Results_Sep28/Orthogroups'
rows=table(op/'Orthogroups.GeneCount.tsv'); cols=[k for k in rows[0] if k not in ['Orthogroup','Total']]
arr=np.array([[int(r[c]) for c in cols] for r in rows]);occ=(arr>0).sum(axis=1)
core_ids=[r['Orthogroup'] for r,ok in zip(rows,(arr==1).all(axis=1)) if ok]
E['orthogroups']={'path':str(op),'samples':len(cols),'columns':cols,'orthogroups':len(rows),
 'strict_single_copy_all':len(core_ids),'present_in_all':int((occ==len(cols)).sum()),
 'present_ge95pct':int((occ>=np.ceil(.95*len(cols))).sum()),'single_sample_orthogroups':int((occ==1).sum()),
 'occupancy_histogram':dict(collections.Counter(map(int,occ))),
 'per_sample':[{'sample':norm(c),'assigned_gene_copies':int(arr[:,i].sum()),'orthogroups_present':int((arr[:,i]>0).sum()),'multicopy_groups':int((arr[:,i]>1).sum())} for i,c in enumerate(cols)],
 'selected_core_file_matches':set((op/'core.txt').read_text().split())==set(core_ids)}
members={r['Orthogroup']:r for r in table(op/'Orthogroups.tsv')}
mafft=ROOT/'results/0_4_mafft_iqtree/Sep28/mafft'
files=sorted(mafft.glob('*_mafft_renamed.fa'))
bad=[]; lengths=[]; ah={}; counts=collections.Counter(); reversed_count=0
for p in files:
    og=p.name.split('_mafft')[0];recs=list(fasta(p)); ls={len(s) for _,s in recs};lengths.append(next(iter(ls)))
    if len(ls)!=1 or len(recs)!=len(cols):bad.append({'file':p.name,'issue':'length/count'})
    seen=set()
    for h,s in recs:
        bits=h.split();c=bits[0]; gene=bits[1] if len(bits)>1 else ''
        if gene.startswith('_R_'):reversed_count+=1;gene=gene[3:]
        expected=members.get(og,{}).get(c,'')
        if gene not in [v.strip() for v in expected.split(',')]:bad.append({'file':p.name,'sample':c,'gene':gene,'expected':expected})
        if c in seen:bad.append({'file':p.name,'issue':'duplicate sample','sample':c})
        seen.add(c);counts[c]+=1;ah.setdefault(c,hashlib.sha256()).update(s.encode())
    if len(lengths)%1000==0:print('Checked alignments',len(lengths),flush=True)
concat=[]
for h,s in fasta(ROOT/bundles['genome']):
    c=h.split()[0]; concat.append({'id':c,'length':len(s),'ambiguous_or_gap_percent':round(100*(len(s)-sum(s.upper().count(b) for b in 'ACGT'))/len(s),3),'matches_lexical_concatenation':c in ah and hashlib.sha256(s.encode()).hexdigest()==ah[c].hexdigest()})
E['concatenation']={'files_checked':len(files),'total_alignment_columns':sum(lengths),'mapping_violation_count':len(bad),'mapping_violation_examples':bad[:30],'reverse_headers':reversed_count,'loci_per_sample':dict(counts),'sequences':concat,'selected_loci_match':set(core_ids)=={p.name.split('_mafft')[0] for p in files}}
print('Completed concatenation audit',flush=True)

# Marker availability is sequence occupancy, not verified orthology.
E['markers']={}
for p in (ROOT/'results/Houbraken_genes/merged/individual').glob('*.fa'):
    rr=list(fasta(p));out=[]
    for h,s in rr:
        valid=sum(s.upper().count(c) for c in 'ACGT')
        out.append({'id':h,'sample':norm(h.split()[0]),'aligned_length':len(s),'acgt_bases':valid,'non_acgt_percent':round(100*(len(s)-valid)/len(s),2)})
    E['markers'][p.stem]={'path':str(p),'records':len(rr),'unique_sample_ids':len(set(r['sample'] for r in out)),'rows':out}

# QC snapshot, data availability and current cohort.
qc=ROOT/'terreus_qc_check';E['qc']={}
for fn in ['complete_genome_qc_summary.tsv','complete_transcript_qc_summary.tsv','data_availability.tsv']:
    rr=table(qc/'complete_data'/fn);E['qc'][fn]={'rows':rr,'count':len(rr)}
gen=set(E['trees']['genome']['normalized_tips']);qg={r['sample'] for r in E['qc']['complete_genome_qc_summary.tsv']['rows']}
E['qc']['cohort_comparison_exact_names']={'tree_not_qc':sorted(gen-qg),'qc_not_tree':sorted(qg-gen),'overlap':len(gen&qg)}
E['metadata']=table(ROOT/'config/samples.tsv')

# Current stage output coverage, based on specific artifacts, not just folders.
patterns={
 'spades':['*/*/scaffolds.fasta'], 'fastp':['*/*.json'],
 'busco':['*/short_summary.specific.*.json'],
 'annotation':['*/*/funannotate/annotate_results/*.summary.json','*/funannotate/annotate_results/*.summary.json'],
 'prediction':['*/*/funannotate/predict_results/*.summary.json','*/funannotate/predict_results/*.summary.json'],
 'signalp':['*/_summary.signalp5'],'dbcan':['*/overview.tsv'],'effectorp':['*/effectorp_res'],
 'variants':['*/*/variants_atLeast3PASS*.vcf'],'sv':['*/*/SV_calling_finished.txt','*/*/*.vcf'],
 'cnv':['*/*/*.vcf'],'integrated':['*/*/SV_and_CNV_variant_calling.vcf'],
 'annotated_sv':['*/*/*.tab'], 'optimized':['*/*/optimized_parameters.json']}
stages={'spades':'0_1_spades','fastp':'0_0_fastp','busco':'0_2_2_busco','annotation':'0_2_repeat_mask_funannotate','prediction':'0_2_repeat_mask_funannotate','signalp':'0_7_signalp','dbcan':'0_8_dbcan','effectorp':'0_9_effectorp','variants':'1_3_perSVade_small_vars','sv':'1_8_perSVade_SVs','cnv':'1_9_perSVade_call_CNVs','integrated':'1_10_perSVade_integrate','annotated_sv':'1_11_perSVade_SV+CNV_annotated','optimized':'1_7_perSVade_optimized_parameters'}
E['coverage']={}
for key,pats in patterns.items():
    fs=sorted(set(p for pat in pats for p in (ROOT/'results'/stages[key]).glob(pat) if p.is_file()))
    E['coverage'][key]={'patterns':pats,'count':len(fs),'files':[str(p.relative_to(ROOT)) for p in fs]}

# Functional output sanity checks; preserve individual values and IDs.
E['signalp']=[]
for p in (ROOT/'results/0_7_signalp').glob('*/_summary.signalp5'):
    ls=p.read_text().splitlines();rr=list(csv.DictReader(ls[1:],delimiter='\t'));motifs=[]
    for row in rr:
        m=re.search(r'\. ([A-Z]+-[A-Z]+)\.',row.get('CS Position',''))
        if m:motifs.append(m.group(1).replace('-',''))
    E['signalp'].append({'sample':p.parent.name,'n':len(rr),'predictions':dict(collections.Counter(x['Prediction'] for x in rr)), 'cleavage_motifs':len(motifs),'motifs_only_ACGT':sum(set(m)<=set('ACGT') for m in motifs),'header':ls[0]})
E['effectorp_id_truncation']=[]
for p in (ROOT/'results/0_9_effectorp').glob('*/effectorp_res'):
    rr=table(p);ids=[x['# Identifier'] for x in rr if x.get('# Identifier')];badids=[x for x in ids if x[:14]!=x.split()[0]]
    E['effectorp_id_truncation'].append({'sample':p.parent.name,'rows':len(ids),'truncated_different_from_first_token':len(badids),'example':badids[:1]})

# Extract existing embedded notebook figures for visual inspection (no execution of notebooks).
images=[];nbinfo=[]
for p in (ROOT/'notebooks').glob('*.ipynb'):
    nb=json.loads(p.read_text());errors=[]
    for i,c in enumerate(nb['cells']):
        for j,o in enumerate(c.get('outputs',[])):
            if o.get('output_type')=='error':errors.append({'cell':i,'ename':o.get('ename'),'evalue':o.get('evalue')})
            if 'image/png' in o.get('data',{}):
                dest=OUT/f'{p.stem}_cell{i}_output{j}.png';v=o['data']['image/png'];dest.write_bytes(base64.b64decode(''.join(v) if isinstance(v,list) else v));images.append(str(dest))
    nbinfo.append({'path':str(p),'cells':len(nb['cells']),'stored_errors':errors})
E['notebooks']=nbinfo;E['extracted_existing_figures']=images

(OUT/'audit_evidence.json').write_text(json.dumps(E,indent=2),encoding='utf-8')
print(json.dumps({'report':str(OUT/'audit_evidence.json'),'og':{k:v for k,v in E['orthogroups'].items() if k not in ['columns','per_sample','occupancy_histogram']},'concatenation':{k:v for k,v in E['concatenation'].items() if k not in ['sequences','loci_per_sample']},'coverage':{k:v['count'] for k,v in E['coverage'].items()},'cohorts':E['qc']['cohort_comparison_exact_names'],'notebook_errors':nbinfo,'images':images},indent=2))
