from pathlib import Path
import json,re,collections,hashlib
BASE=Path(__file__).resolve().parent;OUT=BASE/'taxonomy_reference_build/2026-11-10'
LIB=OUT/'reference_library';WORK=OUT/'analysis';WORK.mkdir(exist_ok=True)
rows=json.loads((LIB/'reference_rows.json').read_text(encoding='utf-8'))
acc=json.loads((LIB/'accepted_species.json').read_text(encoding='utf-8'))
def fasta(p):
    d={};name=None
    for ln in p.read_text().splitlines():
        if ln.startswith('>'):name=ln[1:].split()[0];assert name not in d;d[name]=''
        elif name:d[name]+=ln.strip().upper()
    return d
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()),encoding='utf-8')
conflicts=[]
for d in rows:
    matches=[a for a in acc if all(str(a.get(l+'GB') or '')==str(d.get(l+'GB') or '') for l in ['BenA','CaM'] if a.get(l+'GB'))]
    d['accepted_type_name']=matches[0]['ID'].replace(' ','_') if len(matches)==1 else ''
    d['analysis_species']=d['accepted_type_name'] or d['species']
    d['type_evidence']='accepted_list_accession_match' if d['accepted_type_name'] else 'reference_table_flag' if d['type_flag'] else ''
    if d['accepted_type_name'] and d['accepted_type_name']!=d['species']:conflicts.append({k:d[k] for k in ['reference_id','species','accepted_type_name','type_flag']})
# Keep all Terrei strains plus four representative Flavipedes ex-types.
outgroup_species=['Aspergillus_flavipes','Aspergillus_spelaeus','Aspergillus_neoflavipes','Aspergillus_neoniveus']
selected=[d for d in rows if d['section']=='_sect_Terrei' or (d['type_flag'] and d['species'] in outgroup_species)]
for d in selected:d['role']='ingroup_reference' if d['section']=='_sect_Terrei' else 'outgroup_reference'
(WORK/'reference_metadata.json').write_text(json.dumps(selected,indent=2,ensure_ascii=False),encoding='utf-8')
(LIB/'taxonomy_conflicts.json').write_text(json.dumps(conflicts,indent=2),encoding='utf-8')
missing=[a['ID'] for a in acc if a['SECTION']=='Terrei' and not any(d['accepted_type_name']==a['ID'].replace(' ','_') for d in selected)]
assert not missing,missing
localmeta={};src=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\Houbraken_genes\merged\individual')
for locus,fn in [('BenA','BENA56_aligned.fa'),('CaM','CAM1_aligned.fa'),('RPB2','RPB2_aligned.fa')]:
    local=fasta(src/fn);clean={}
    for name,seq in local.items():
        sample=re.sub(r'^(?:terreus|hortae|pseudoterreus|alabamensis|aureoterreus|floccosus)_','',name).removesuffix('_transcripts')
        ident='LOCAL_'+sample
        m=localmeta.setdefault(ident,{'sample':sample,'original_labels':[],'loci':{},'role':'project_sample'})
        if name not in m['original_labels']:m['original_labels'].append(name)
        seq=seq.replace('-','').replace('?','N')
        assert not set(seq)-set('ACGTRYKMSWBDHVN'),(name,locus)
        m['loci'][locus]={'file':str(src/fn),'ungapped_length':len(seq),'ACGT':sum(seq.count(c) for c in 'ACGT'),'file_sha256':hashlib.sha256((src/fn).read_bytes()).hexdigest()}
        if any(c in seq for c in 'ACGT'):clean[ident]=seq
    writefa(WORK/(locus+'_local_ungapped.fasta'),clean)
    writefa(WORK/(locus+'_references.fasta'),{d['reference_id']:d[locus+'_sequence'] for d in selected if d[locus+'_sequence']})
    writefa(WORK/(locus+'_type_queries.fasta'),{d['reference_id']:d[locus+'_sequence'] for d in selected if d['accepted_type_name'] and d[locus+'_sequence']})
    print(locus,'local',len(clean),'refs',sum(bool(d[locus+'_sequence']) for d in selected))
(WORK/'local_metadata.json').write_text(json.dumps(localmeta,indent=2),encoding='utf-8')
print('LOCAL',len(localmeta),'REFERENCES',len(selected),'TYPE_NAME_CONFLICTS',conflicts)
