"""Download linked publications, taxonomic records and full source metadata."""
from harvest_ncbi_core_metadata import *
import collections

def load(name):return json.loads((OUT/name).read_text(encoding='utf-8'))
indices=load('entrez_indices.json');rows=load('tip_source_mapping.json')
assemblies=load('assembly_reports.json');packages=load('sra_packages_parsed.json')
links=[];errors=[]
for db in ['assembly','sra','biosample','bioproject']:
    try:
        root=eu('elink',{'dbfrom':db,'db':'pubmed','id':indices[db]['uids'],'retmode':'xml'},db+'_pubmed_links.xml','xml')
        for ls in root.findall('LinkSet'):
            src=[x.text for x in ls.findall('IdList/Id')]
            for dest in ls.findall('LinkSetDb'):
                for uid in dest.findall('Link/Id'):
                    links.append({'source_db':db,'source_uids':src,'link_name':text(dest,'LinkName'),'pmid':uid.text,'evidence':'NCBI Entrez cross-link'})
        print('LINK',db,len(links),flush=True)
    except Exception as exc:errors.append({'stage':'pubmed_elink','db':db,'error':str(exc)})

bio=ET.parse(RAW/'biosample_all.xml').getroot()
biosamples=[]
for el in bio.findall('.//BioSample'):
    biosamples.append({'accession':el.get('accession'),'uid':el.get('id'),
      'title':text(el,'Description/Title'),'organism':el.find('Description/Organism').attrib if el.find('Description/Organism') is not None else {},
      'attributes':[dict(a.attrib,value=''.join(a.itertext()).strip()) for a in el.findall('Attributes/Attribute')],
      'ids':[dict(a.attrib,value=''.join(a.itertext()).strip()) for a in el.findall('Ids/Id')],
      'dates':{k:v for k,v in el.attrib.items() if k not in ['id','accession']},
      'owner':text(el,'Owner/Name'),'links':[dict(a.attrib,value=''.join(a.itertext()).strip()) for a in el.findall('.//Link')]})
save('biosample_records.json',biosamples)
bp=ET.parse(RAW/'bioproject_all.xml').getroot()
projects=[]
for el in bp.findall('.//DocumentSummary'):
    archive=el.find('.//ProjectID/ArchiveID')
    if archive is None:continue
    projects.append({'accession':archive.get('accession'),'uid':el.get('uid') or archive.get('id'),
      'title':text(el,'.//ProjectDescr/Title'),'description':text(el,'.//ProjectDescr/Description'),
      'publications':[dict(x.attrib,db_type=text(x,'DbType'),title=text(x,'StructuredCitation/Title'),text=''.join(x.itertext()).strip()) for x in el.findall('.//Publication')],
      'external_links':[xmlstr(x) for x in el.findall('.//ExternalLink')]})
    for x in el.findall('.//Publication'):
        if (x.get('dbType') or text(x,'DbType'))=='ePubmed' and x.get('id'):
            links.append({'source_db':'bioproject','source_uids':[el.get('uid') or archive.get('id')],
                          'source_accession':archive.get('accession'),'pmid':x.get('id'),'evidence':'BioProject publication field'})
save('bioproject_records.json',projects)
for r in rows:
    for pmid in r.get('local_genbank_pubmed',[]):links.append({'source_db':'local_genbank','source_sample':r['sample_id'],'source_accessions':r['assembly_accessions'],'pmid':pmid,'evidence':'Local downloaded GenBank REFERENCE/PUBMED field'})
for p in packages:
    for l in p['study_links']:
        if l['db'].lower()=='pubmed' and l['id'].isdigit():links.append({'source_db':'sra_study','source_runs':p['run_accessions'],'pmid':l['id'],'evidence':'SRA study reference'})

# Search exact BioProject identifiers in PubMed to catch publications absent from cross-links.
# Keep these as candidate associations, distinct from explicit database links.
for pr in projects:
    if pr['accession'] not in indices['bioproject']['requested_accessions']:continue
    if any(l.get('source_accession')==pr['accession'] or (l['source_db']=='bioproject' and str(pr['uid']) in l.get('source_uids',[])) for l in links):continue
    try:
        found=search('pubmed',[pr['accession']],'pubmed_search_'+pr['accession'])
        for pmid in found:links.append({'source_db':'bioproject','source_uids':[pr['uid']],'source_accession':pr['accession'],'pmid':pmid,'evidence':'PubMed exact-accession search; association requires record review'})
    except Exception as exc:errors.append({'stage':'pubmed_search','accession':pr['accession'],'error':str(exc)})

save('publication_links.json',links)
pmids=sorted(set(l['pmid'] for l in links))
if pmids:
    pub=eu('efetch',{'db':'pubmed','id':','.join(pmids),'retmode':'xml'},'pubmed_publications.xml','xml')
    citations=[]
    for art in pub.findall('PubmedArticle'):
        cites=art.find('MedlineCitation');a=cites.find('Article')
        citations.append({'pmid':text(cites,'PMID'),'title':text(a,'ArticleTitle'),'journal':text(a,'Journal/Title'),
          'year':text(a,'Journal/JournalIssue/PubDate/Year') or text(a,'Journal/JournalIssue/PubDate/MedlineDate'),
          'authors':[' '.join(filter(None,[text(x,'LastName'),text(x,'Initials')])) or text(x,'CollectiveName') for x in a.findall('AuthorList/Author')],
          'doi':next((x.text for x in art.findall('.//ArticleId') if x.get('IdType')=='doi'),''),
          'pmc':next((x.text for x in art.findall('.//ArticleId') if x.get('IdType')=='pmc'),''),
          'url':'https://pubmed.ncbi.nlm.nih.gov/'+text(cites,'PMID')+'/',
          'abstract':[{'label':x.get('Label',''),'text':''.join(x.itertext())} for x in a.findall('Abstract/AbstractText')]})
    save('publication_records.json',citations)
    print('PUBLICATIONS',len(citations),flush=True)

# Taxonomy is downloaded as reported metadata, not an independent taxonomic adjudication.
taxids=sorted(set(str(b['organism'].get('taxonomy_id','')) for b in biosamples)-{''})
if taxids:eu('efetch',{'db':'taxonomy','id':','.join(taxids),'retmode':'xml'},'taxonomy_records.xml','xml')

# Record-level assembly cross-links can reveal deposited versions of the same BioSample.
# These are related records, never automatically substituted for the assembly used in the tree.
try:
    ar=eu('elink',{'dbfrom':'biosample','db':'assembly','id':indices['biosample']['uids'],'retmode':'xml'},'biosample_assembly_links.xml','xml')
    related=[]
    for ls in ar.findall('LinkSet'):
        related.append({'biosample_uids':[x.text for x in ls.findall('IdList/Id')],
                        'assembly_uids':sorted(set(x.text for x in ls.findall('LinkSetDb/Link/Id')))})
    save('related_assembly_links.json',related)
    uids=sorted({a for r in related for a in r['assembly_uids']})
    if uids:eu('esummary',{'db':'assembly','id':','.join(uids),'retmode':'json'},'related_assemblies_summary.json')
    print('RELATED ASSEMBLIES',len(uids),flush=True)
except Exception as exc:errors.append({'stage':'related_assemblies','error':str(exc)})

# Local-only identifiers are not sent to external services. Their NCBI linkage
# remains unresolved unless the user supplies an accession or authorizes a query.
save('enrichment_errors.json',errors)
print('ENRICHMENT DONE',len(errors),'errors',flush=True)
