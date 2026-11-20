from terreus_reference_sources import *
import csv
L=OUT/'reference_library'
rows=list(csv.DictReader((L/'marker_associated_publications.tsv').open(),delimiter='\t'));pmids=sorted({r['PMID'] for r in rows if r['PMID']})
root=eu('efetch',{'db':'pubmed','id':','.join(pmids),'retmode':'xml'},'marker_publications_pubmed.xml')
pubs=[]
for a in root.findall('.//PubmedArticle'):
    art=a.find('./MedlineCitation/Article');title=''.join(art.find('ArticleTitle').itertext());pmid=a.findtext('./MedlineCitation/PMID');doi=next((x.text for x in a.findall('./PubmedData/ArticleIdList/ArticleId') if x.get('IdType')=='doi'),'')
    authors=[' '.join(filter(None,[x.findtext('LastName'),x.findtext('Initials')])) or x.findtext('CollectiveName','') for x in art.findall('./AuthorList/Author')]
    pubs.append({'PMID':pmid,'title':title,'authors':authors,'journal':art.findtext('./Journal/Title'),'year':art.findtext('./Journal/JournalIssue/PubDate/Year') or art.findtext('./Journal/JournalIssue/PubDate/MedlineDate'),'doi':doi,'url':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/'})
(L/'marker_publications_enriched.json').write_text(json.dumps(pubs,indent=2,ensure_ascii=False),encoding='utf-8')
with (L/'marker_publications.ris').open('w',encoding='utf-8') as f:
    for d in pubs:
        f.write('TY  - JOUR\nTI  - '+d['title']+'\n')
        for au in d['authors']:f.write('AU  - '+au+'\n')
        for tag,key in [('JO','journal'),('PY','year'),('DO','doi'),('UR','url')]:
            if d[key]:f.write(tag+'  - '+d[key]+'\n')
        f.write('ER  - \n\n')
print('PUBLICATIONS ENRICHED',len(pubs))
