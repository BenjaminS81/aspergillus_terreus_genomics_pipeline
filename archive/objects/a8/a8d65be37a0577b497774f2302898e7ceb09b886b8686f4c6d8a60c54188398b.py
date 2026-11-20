from pathlib import Path
import openpyxl,json,re,collections,csv
BASE=Path(__file__).resolve().parent
OUT=BASE/'taxonomy_reference_build/2026-11-10'
SRC=OUT/'sources'
LIB=OUT/'reference_library';LIB.mkdir(exist_ok=True)
accepted=[]
w=openpyxl.load_workbook(SRC/'ICPA_accepted_species_20261001.xlsx',read_only=True,data_only=True)
s=w.worksheets[1];rows=s.iter_rows(values_only=True);head=next(rows)
for n,row in enumerate(rows,2):
    d=dict(zip(head,row))
    if d.get('SECTION') in ['Terrei','Flavipedes'] or d.get('ID') in ['Aspergillus candidus','Aspergillus fumigatus']:
        d['source_row']=n;accepted.append(d)
(LIB/'accepted_species.json').write_text(json.dumps(accepted,indent=2,ensure_ascii=False),encoding='utf-8')
ref=[]
w=openpyxl.load_workbook(SRC/'ICPA_reference_sequences_20261001.xlsx',read_only=True,data_only=True)
s=w['Sheet1'];rows=s.iter_rows(max_col=21,values_only=True);head=next(rows)
for n,row in enumerate(rows,2):
    d=dict(zip(head,row))
    if d.get('section') not in ['_sect_Terrei','_sect_Flavipedes'] and d.get('species') not in ['Aspergillus_candidus','Aspergillus_fumigatus']:continue
    d['source_row']=n
    d['reference_id']='REF_'+str(n)+'_'+re.sub(r'[^A-Za-z0-9_-]','',str(d['species']))+'_'+re.sub(r'[^A-Za-z0-9_-]','',str(d['strain']))
    d['type_flag']=d.get('_2')=='T_'
    for locus in ['BenA','CaM','RPB2']:
        seq=str(d.get(locus+'_Seq_trimmed') or '').strip().upper()
        if seq.startswith('*'):seq=seq[1:]
        if seq in ['','N/A','NA','NONE','-']:seq=''
        if seq and set(seq)-set('ACGTURYKMSWBDHVN-?'):raise ValueError((n,locus,seq[:100],set(seq)-set('ACGTURYKMSWBDHVN-?')))
        d[locus+'_sequence']=seq.replace('-','').replace('?','N').replace('U','T')
    ref.append(d)
(LIB/'reference_rows.json').write_text(json.dumps(ref,indent=2,ensure_ascii=False),encoding='utf-8')
for locus in ['BenA','CaM','RPB2']:
    with (LIB/(locus+'_all_references.fasta')).open('w') as f:
        for d in ref:
            if d[locus+'_sequence']:f.write('>'+d['reference_id']+'\n'+d[locus+'_sequence']+'\n')
cols=['reference_id','species','strain','type_flag','subgenus','section','series','source_row','BenAGB','CaMGB','RPB2GB']
with (LIB/'reference_inventory.tsv').open('w',newline='',encoding='utf-8') as f:
    wr=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore');wr.writeheader();wr.writerows(ref)
print('ROWS',len(ref),'SECTIONS',collections.Counter(d['section'] for d in ref))
for sec in ['_sect_Terrei','_sect_Flavipedes']:
    for d in ref:
        if d['section']==sec:print(d['source_row'],d['species'],d['strain'],d['type_flag'],*[str(d.get(l+'GB'))+':'+str(len(d[l+'_sequence'])) for l in ['BenA','CaM','RPB2']])
print('DISTANT',[(d['reference_id'],d['type_flag']) for d in ref if d['section'] not in ['_sect_Terrei','_sect_Flavipedes']])
