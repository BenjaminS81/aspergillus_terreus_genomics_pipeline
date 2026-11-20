from harvest_ncbi_core_metadata import *
indices=json.loads((OUT/'entrez_indices.json').read_text(encoding='utf-8'))
r=eu('elink',{'dbfrom':'biosample','db':'sra','id':indices['biosample']['uids'],'retmode':'xml'},'biosample_sra_links.xml','xml')
links=[]
for ls in r.findall('LinkSet'):
    links.append({'biosample_uids':[x.text for x in ls.findall('IdList/Id')],
                  'sra_uids':sorted(set(x.text for x in ls.findall('LinkSetDb/Link/Id')))})
save('related_sra_links.json',links)
uids=sorted({u for l in links for u in l['sra_uids']})
print('RELATED SRA EXPERIMENTS',len(uids),flush=True)
assert len(uids)<=1000,'Unexpectedly broad related-record set; inspect before fetching.'
if uids:
    for i in range(0,len(uids),100):
        batch=uids[i:i+100]
        eu('efetch',{'db':'sra','id':','.join(batch),'retmode':'xml'},f'related_sra_{i//100}.xml','xml')
        eu('efetch',{'db':'sra','id':','.join(batch),'rettype':'runinfo','retmode':'text'},f'related_sra_runinfo_{i//100}.csv','text')
summ=json.loads((RAW/'related_assemblies_summary.json').read_text(encoding='utf-8'))['result']
accs=sorted({summ[u]['assemblyaccession'] for u in summ['uids']})
if accs:
    get('https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+','.join(accs)+'/dataset_report?page_size=1000','related_assemblies_dataset_report.json')
eu('efetch',{'db':'pmc','id':'11278011','retmode':'xml'},'section_terrei_publication.xml','xml')
print('RELATED RECORDS DONE',flush=True)
