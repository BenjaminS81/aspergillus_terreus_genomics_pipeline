from harvest_ncbi_core_metadata import *

def load(n):return json.loads((OUT/n).read_text(encoding='utf-8'))
indices=load('entrez_indices.json');rows=load('tip_source_mapping.json');pkgs=load('sra_packages_parsed.json');ars=load('assembly_reports.json')
eu('esummary',{'db':'sra','id':','.join(indices['sra']['uids']),'retmode':'json'},'sra_tree_runs_summary.json')
summary=load('raw/assembly_all_summary.json')['result']
errors=[]
for uid in summary['uids']:
    a=summary[uid];acc=a['assemblyaccession']
    for key,kind in [('ftppath_assembly_rpt','assembly_report'),('ftppath_stats_rpt','assembly_stats')]:
        url=a.get(key)
        if not url:continue
        try:get(url.replace('ftp://','https://'),acc+'_'+kind+'.txt','text')
        except Exception as exc:errors.append({'stage':kind,'accession':acc,'error':str(exc)})
    print('ASSEMBLY REPORTS',acc,flush=True)

runpkg={run:p for p in pkgs for run in p['run_accessions']}
amap={a['accession']:a for a in ars}
queries={}
for row in rows:
    terms=set(row['assembly_accessions']+row['run_accessions'])
    for run in row['run_accessions']:terms.update(runpkg[run]['bioprojects'])
    for acc in row['assembly_accessions']:
        bp=amap.get(acc,{}).get('assembly_info',{}).get('bioproject_accession')
        if bp:terms.add(bp)
    if terms:queries[tuple(sorted(terms))]=queries.get(tuple(sorted(terms)),[])+[row['sample_id']]
hits=[]
for i,(terms,samples) in enumerate(queries.items(),1):
    label='pmc_search_'+hashlib.sha256('|'.join(terms).encode()).hexdigest()[:16]
    try:
        ids=search('pmc',list(terms),label)
        hits.append({'samples':samples,'public_accession_terms':list(terms),'pmc_uids':ids,'evidence':'PMC exact public-accession search; candidate mention, not yet an original-study attribution'})
        if ids:print('PMC HITS',samples,ids,flush=True)
    except Exception as exc:errors.append({'stage':'pmc_search','samples':samples,'error':str(exc)})
save('pmc_accession_search_hits.json',hits)
uids=sorted({x for h in hits for x in h['pmc_uids']})
if uids:
    eu('esummary',{'db':'pmc','id':','.join(uids),'retmode':'json'},'pmc_candidates_summary.json')
    eu('elink',{'dbfrom':'pmc','db':'pubmed','id':uids,'retmode':'xml'},'pmc_candidates_pubmed_links.xml','xml')

# Targeted title searches identified by public-accession web searches.
titles=['Whole-genome sequencing phenotypic characterization antifungal susceptibility profiles three Aspergillus hortae clinical isolates Colombia',
        'Species Identification Orthologous Allergen Prediction Expression Genus Aspergillus',
        'Genomic Analysis Aspergillus Section Terrei Reveals High Potential Secondary Metabolite Production Plant Biomass Degradation']
titlehits=[]
for i,title in enumerate(titles):
    d=eu('esearch',{'db':'pubmed','term':title,'retmode':'json','retmax':50},'target_title_'+str(i)+'.json')
    titlehits.append({'query':title,'pmids':d['esearchresult']['idlist']})
save('targeted_title_searches.json',titlehits)
pubids=sorted({x for t in titlehits for x in t['pmids']})
if uids:
    links=ET.parse(RAW/'pmc_candidates_pubmed_links.xml').getroot()
    pubids=sorted(set(pubids)|{x.text for block in links.findall('.//LinkSetDb') if text(block,'LinkName')=='pmc_pubmed' for x in block.findall('Link/Id')})
if pubids:eu('efetch',{'db':'pubmed','id':','.join(pubids),'retmode':'xml'},'supplementary_pubmed_publications.xml','xml')
save('supplementary_errors.json',errors)
print('SUPPLEMENT DONE',len(uids),'PMC candidates;',len(pubids),'PubMed candidates;',len(errors),'errors',flush=True)
