from harvest_ncbi_core_metadata import *
queries={
 'M6925':'"Whole-Genome Sequencing of Aspergillus terreus Species Complex"[Title]',
 'TT0037':'Aspergillus AND (TT0037 OR Hypotaenidia OR Okinawa)',
 'w25':'Aspergillus AND "w25"',
 'Brazil_study':'Aspergillus AND terreus AND ("PRJNA554775" OR "SRR9691930" OR "SRR9691933" OR "SRR9691934" OR "SRR9691935" OR "SRR9691936")'}
hits=[]
for label,q in queries.items():
    d=eu('esearch',{'db':'pubmed','term':q,'retmode':'json','retmax':100},'target_public_'+label+'.json')
    hits.append({'label':label,'query':q,'pmids':d['esearchresult']['idlist']})
save('targeted_public_strain_searches.json',hits)
pmids=sorted({p for h in hits for p in h['pmids']})
if pmids:eu('efetch',{'db':'pubmed','id':','.join(pmids),'retmode':'xml'},'targeted_public_strain_pubmed.xml','xml')
print(hits,flush=True)
