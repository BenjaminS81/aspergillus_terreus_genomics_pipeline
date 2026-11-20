from terreus_reference_sources import *
import re,collections,csv
rows=json.loads((OUT/'reference_library/reference_rows.json').read_text(encoding='utf-8'))
accs=sorted({str(d.get(l+'GB')) for d in rows for l in ['BenA','CaM','RPB2'] if re.fullmatch(r'[A-Z]{1,2}\d{5,8}(?:\.\d+)?',str(d.get(l+'GB')))})
records={}
for i in range(0,len(accs),40):
    batch=accs[i:i+40]
    r=eu('efetch',{'db':'nuccore','id':','.join(batch),'rettype':'gb','retmode':'xml'},f'genbank_markers_{i:04d}.xml')
    for gb in r.findall('.//GBSeq'):
        d={k:gb.findtext('GBSeq_'+k) for k in ['primary-accession','accession-version','definition','organism','sequence','length','taxonomy']}
        d['references']=[{k:x.findtext('GBReference_'+k) for k in ['title','journal','pubmed']} for x in gb.findall('./GBSeq_references/GBReference')]
        d['source_qualifiers']={q.findtext('GBQualifier_name'):q.findtext('GBQualifier_value') for f in gb.findall('./GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key')=='source' for q in f.findall('./GBFeature_quals/GBQualifier')}
        records[d['primary-accession']]=d
    print('FETCH',min(i+40,len(accs)),'/',len(accs),'records',len(records),flush=True)
save('reference_library/genbank_records.json',records)
comp=str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN')
validation=[]
for d in rows:
    for l in ['BenA','CaM','RPB2']:
        a=str(d.get(l+'GB') or '');seq=d[l+'_sequence'];gb=records.get(a.split('.')[0]);status='no_sequence'
        if seq:
            if a.startswith('GCA'):status='assembly_derived_curated_sequence_not_nuccore_accession'
            elif gb:
                full=(gb['sequence'] or '').upper()
                status='exact_subsequence' if seq in full else 'reverse_complement_subsequence' if seq.translate(comp)[::-1] in full else 'sequence_discrepancy_review'
            else:status='accession_not_retrieved'
        validation.append({'reference_id':d['reference_id'],'locus':l,'accession':a,'version':gb['accession-version'] if gb else '', 'sequence_length':len(seq),'status':status,'genbank_organism':gb['organism'] if gb else ''})
save('reference_library/genbank_validation.json',validation)
with (OUT/'reference_library/genbank_validation.tsv').open('w',newline='',encoding='utf-8') as f:
    wr=csv.DictWriter(f,fieldnames=list(validation[0]),delimiter='\t');wr.writeheader();wr.writerows(validation)
print('VALIDATION',dict(collections.Counter(x['status'] for x in validation)),flush=True)
