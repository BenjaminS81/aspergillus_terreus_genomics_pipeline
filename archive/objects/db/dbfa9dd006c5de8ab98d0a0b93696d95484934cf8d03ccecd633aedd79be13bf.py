from terreus_reference_sources import *
import collections
accessions=[f'FJ{i}' for i in range(491703,491732)]+[f'FJ{i}' for i in range(531192,531244)]
out=OUT/'reference_library/historical_2011';out.mkdir(exist_ok=True)
records=[]
for i in range(0,len(accessions),40):
    r=eu('efetch',{'db':'nuccore','id':','.join(accessions[i:i+40]),'rettype':'gb','retmode':'xml'},f'Samson_2011_accession_range_{i:03d}.xml')
    for gb in r.findall('.//GBSeq'):
        d={k:gb.findtext('GBSeq_'+k) for k in ['primary-accession','accession-version','definition','organism','sequence','length']}
        d['source_qualifiers']={q.findtext('GBQualifier_name'):q.findtext('GBQualifier_value') for f in gb.findall('./GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key')=='source' for q in f.findall('./GBFeature_quals/GBQualifier')}
        txt=d['definition'].lower();d['locus']='BenA' if 'tubulin' in txt else 'CaM' if 'calmodulin' in txt else 'ITS' if 'internal transcribed spacer' in txt else 'other';records.append(d)
(out/'records.json').write_text(json.dumps(records,indent=2))
for locus in ['BenA','CaM','ITS','other']:
    (out/(locus+'.fasta')).write_text(''.join('>'+d['accession-version']+'\n'+d['sequence'].upper()+'\n' for d in records if d['locus']==locus))
print('HISTORICAL 2011',len(records),collections.Counter(d['locus'] for d in records))
