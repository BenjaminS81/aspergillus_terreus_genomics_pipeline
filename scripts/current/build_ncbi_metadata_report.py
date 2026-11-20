"""Build a locally joined metadata catalogue and a conservative publication map."""
from pathlib import Path
import json,csv,re,html,collections,xml.etree.ElementTree as ET,hashlib,zipfile
BASE=Path(__file__).resolve().parent
OUT=BASE/'review_reports/2026-11-10/ncbi_core_tree_metadata';RAW=OUT/'raw'
def load(n):return json.loads((OUT/n).read_text(encoding='utf-8'))
def save(n,obj):(OUT/n).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
def tx(el,path):
    n=el.find(path);return ''.join(n.itertext()).strip() if n is not None else ''
def tsv(n,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/n).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def joined(v):return '; '.join(str(x) for x in v)
rows=load('tip_source_mapping.json');ars={a['accession']:a for a in load('assembly_reports.json')}
bs={b['accession']:b for b in load('biosample_records.json')}
bp={p['accession']:p for p in load('bioproject_records.json')}
pk={r:p for p in load('sra_packages_parsed.json') for r in p['run_accessions']}
runinfo={r['Run']:r for r in load('runinfo_records.json')}
asu=load('raw/assembly_all_summary.json')['result']
asmuid={u:set(filter(None,[asu[u].get('assemblyaccession'),asu[u].get('synonym',{}).get('genbank'),asu[u].get('synonym',{}).get('refseq')])) for u in asu['uids']}
sras=load('raw/sra_tree_runs_summary.json')['result']
srauid={u:re.findall(r'acc="([SED]RR\d+)"',sras[u]['runs']) for u in sras['uids']}
buid={b['uid']:b['accession'] for b in bs.values()};puid={str(p['uid']):p['accession'] for p in bp.values()}

publications={}
for fname in ['pubmed_publications.xml','supplementary_pubmed_publications.xml','reference_followup_pubmed.xml','verified_followup_pubmed.xml']:
    for art in ET.parse(RAW/fname).getroot().findall('PubmedArticle'):
        a=art.find('MedlineCitation/Article');pmid=tx(art,'MedlineCitation/PMID')
        publications[pmid]={'pmid':pmid,'title':tx(a,'ArticleTitle'),'journal':tx(a,'Journal/Title'),
          'year':tx(a,'Journal/JournalIssue/PubDate/Year') or tx(a,'Journal/JournalIssue/PubDate/MedlineDate'),
          'authors':[' '.join(filter(None,[tx(x,'LastName'),tx(x,'Initials')])) or tx(x,'CollectiveName') for x in a.findall('AuthorList/Author')],
          'doi':next((x.text for x in art.findall('.//ArticleId') if x.get('IdType')=='doi'),''),
          'pmc':next((x.text for x in art.findall('.//ArticleId') if x.get('IdType')=='pmc'),''),
          'url':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/',
          'abstract':[{'label':x.get('Label',''),'text':''.join(x.itertext())} for x in a.findall('Abstract/AbstractText')]}

associations=[]
def assoc(pmid,samples,role,evidence):
    assert pmid in publications,pmid
    for sample in samples:
        assert sample in {r['sample_id'] for r in rows},sample
        associations.append({'sample_id':sample,'pmid':pmid,'role':role,'evidence':evidence,'url':publications[pmid]['url']})
assoc('27284150',['45A'],'genome_publication','Assembly/BioProject cross-link; strain also called KM017963.')
assoc('33515286',['ATCC20542'],'genome_publication','Assembly/BioProject link and explicit PRJNA622971 data availability.')
assoc('36604572',['ATCC20541'],'genome_publication','Assembly/BioProject link and explicit PRJNA861866 data availability.')
assoc('31806745',['TN-484','TN-484_DRR191518'],'genome_publication','Original TN-484 genome announcement; source run DRR191518. Local reference is assembly version .2, so do not assume every published statistic describes that version.')
assoc('32299870',['IFO6365','IFO6365_DRR204676'],'genome_publication','Original IFO6365 genome announcement and shared deposited source identifiers.')
assoc('32108289',['M6925'],'genome_publication','Publication explicitly describes genome sequencing of clinical strain M6925; DOI 10.1007/s11046-020-00433-0.')
assoc('36294615',['SRR19634133'],'genome_publication','M7 strain metadata and PRJNA848651 explicitly identified in publication.')
assoc('41055369',['UA-600B','UA-600B_SRR33476766'],'genome_publication','Article table/data availability identifies GCA_049862995.1 and SRR33476766.')
assoc('41701753',['SRR24718917','SRR24718918','SRR24718919'],'genome_and_phenotype_publication','2026 article explicitly identifies PRJNA975750 and BioSamples SAMN35344990/991/992, matching MCA7/8/10.')
assoc('39997392',['SRR24718917','SRR24718918','SRR24718919'],'associated_project_publication','Shared PRJNA975750; 81-isolate identification/allergen study. Check isolate-specific supplementary tables before importing any phenotype.')
assoc('39057392',['SRR10507650','SRR8399219','SRR8399223'],'genome_publication','Section Terrei article explicitly names CBS116.37, IBT12702 and CBS503.65 as sequenced strains; floccosus also has a direct NCBI publication link.')
assoc('42437896',['AT.1'],'genome_publication','Data availability explicitly names GCA_057814615.1, PRJNA1429839 and SAMN56256256.')
assoc('41497419',['TT0037'],'genome_and_phenotype_publication','Article names TT0037 and deposits DRR582540 and BAAFRR010000001–12; matches the reference BioSample/assembly.')
assoc('41704769',['PA2902','PB4404','P3406A','S1101B','PG6503','R6201Q'],'associated_project_publication','PRJNA833221 is linked to this gut-fungal study. The retrieved main text identifies the project; strain-specific inclusion/clinical cohort assignment needs supplementary confirmation.')
assoc('34485249',['ASM-1','ML-44'],'strain_background','2021 study explicitly identifies ASM-1 as a chemically mutagenized derivative of ML-44. This is not established here as the genome-deposition publication.')
assoc('30394102',['ML-44'],'strain_background','Named-strain natural-products study, not a validated citation for the deposited assembly version.')
assoc('38445873',['IFO6365','TN-484'],'taxonomic_reassessment','Table lists GCA_009932835.1 and GCA_009014675.2 under updated species Aspergillus pseudoterreus, versus Aspergillus terreus on NCBI.')
assoc('38445873',['IFO6365_DRR204676','TN-484_DRR191518'],'taxonomic_reassessment_related_source','The named deposited assemblies were reassessed; applicability to these independently assembled run-derived tips follows their shared BioSample and still needs sequence/annotation validation.')
assoc('33515286',['NIH2624','45A','w25','M6925','IFO6365'],'comparative_reuse','Exact assembly accession appears in comparison table; this is reuse, not the original publication for these samples.')
assoc('36604572',['NIH2624','45A','w25','M6925','ASM-1','ML-44','IFO6365','TN-484'],'comparative_reuse','Named comparative-genome panel in the ATCC20541 article; preserve assembly-version differences.')
assoc('38445873',['TN-484_DRR191518','IFO6365_DRR204676'],'source_identity_caveat','Taxonomic label should not be propagated solely by strain name; compare deposited and run-derived assemblies before a release.')

doi_only=[{'doi':'10.1016/0922-338X(95)91272-7','title':'Breeding of Aspergillus terreus mutant TN-484 for itaconic acid production with high yield','year':'1995','samples':['TN-484','TN-484_DRR191518'],'role':'strain_background','evidence':'Explicit DOI citation in BioProject PRJDB8764; historical strain breeding, not genome sequencing.','url':'https://doi.org/10.1016/0922-338X(95)91272-7'}]

# Sources provide unmodified strings; these are not adjudicated clinical or taxonomic labels.
fields=['strain','isolate','sample_name','geo_loc_name','collection_date','isolation_source','host','host_disease','lat_lon','ploidy','sample_type','type_material']
def key(a):return (a.get('harmonized_name') or a.get('attribute_name','')).lower().replace('-','_').replace(' ','_')
wide=[];long=[];full=[]
related=load('raw/related_assemblies_summary.json')['result']
related_by_bs=collections.defaultdict(set)
for uid in related['uids']:
    a=related[uid]
    if a.get('biosampleaccn'):related_by_bs[a['biosampleaccn']].add(a['assemblyaccession'])
related_runs_by_bs=collections.defaultdict(set)
for p in RAW.glob('related_sra_runinfo_*.csv'):
    for row in csv.DictReader(p.open(encoding='utf-8-sig')):
        if row.get('BioSample'):related_runs_by_bs[row['BioSample']].add(row['Run'])
for row in rows:
    runs=row['run_accessions'];accs=row['assembly_accessions'];bi=set();projects=set()
    for run in runs:
        bi.update(pk[run]['biosamples']);projects.update(pk[run]['bioprojects'])
    for acc in accs:
        ai=ars[acc].get('assembly_info',{})
        if ai.get('biosample',{}).get('accession'):bi.add(ai['biosample']['accession'])
        if ai.get('bioproject_accession'):projects.add(ai['bioproject_accession'])
    attrs=collections.defaultdict(list)
    for b in sorted(bi):
        for a in bs[b]['attributes']:
            k=key(a);attrs[k].append(a['value'])
            long.append({'sample_id':row['sample_id'],'biosample':b,'field':k,'submitted_attribute_name':a.get('attribute_name',''),'value':a['value'],'source_url':'https://www.ncbi.nlm.nih.gov/biosample/'+b})
    pubs=[a for a in associations if a['sample_id']==row['sample_id']]
    r={'sample_id':row['sample_id'],'tree_tip':row['tree_tip'],'mapping_status':row['mapping_status'],
       'used_run_accessions':joined(runs),'local_reference_assembly_accessions':joined(accs),'biosample_accessions':joined(sorted(bi)),
       'bioproject_accessions':joined(sorted(projects)),
       'ncbi_organism':joined(sorted({bs[b]['organism'].get('taxonomy_name','') for b in bi})),
       'ncbi_taxids':joined(sorted({bs[b]['organism'].get('taxonomy_id','') for b in bi})),
       **{f:joined(dict.fromkeys(attrs.get(f,[]))) for f in fields},
       'library_strategy':joined(sorted({runinfo[r].get('LibraryStrategy','') for r in runs})),
       'library_layout':joined(sorted({runinfo[r].get('LibraryLayout','') for r in runs})),
       'sequencing_platform':joined(sorted({runinfo[r].get('Platform','') for r in runs})),
       'instrument':joined(sorted({runinfo[r].get('Model','') for r in runs})),
       'related_public_runs':joined(sorted({x for b in bi for x in related_runs_by_bs[b]})),
       'related_public_assemblies':joined(sorted({x for b in bi for x in related_by_bs[b]})),
       'publication_pmids':joined(sorted({p['pmid'] for p in pubs})),
       'publication_roles':joined(sorted({p['role'] for p in pubs})),
       'taxonomy_review_flag':'Published A. pseudoterreus reassessment conflicts with NCBI A. terreus' if row['sample_id'] in ['IFO6365','IFO6365_DRR204676','TN-484','TN-484_DRR191518'] else '',
       'biological_relationship_flag':'ASM-1 is a reported mutagenized derivative of ML-44; not independent environmental sampling' if row['sample_id'] in ['ASM-1','ML-44'] else '',
       'mapping_caveat':row['mapping_caveat']}
    wide.append(r)
    full.append({'summary':r,'local_source_evidence':row,'ncbi_assembly_reports':[ars[a] for a in accs],
                 'ncbi_biosamples':[bs[b] for b in sorted(bi)],'ncbi_bioprojects':[bp[p] for p in sorted(projects) if p in bp],
                 'ncbi_sra_packages':[pk[r] for r in runs],'ncbi_runinfo':[runinfo[r] for r in runs],
                 'publication_associations':pubs})
assert len(wide)==66 and len({r['sample_id'] for r in wide})==66
assert sum(bool(r['biosample_accessions']) for r in wide)==51
assert all(len(r['biosample_accessions'].split('; '))==1 for r in wide if r['biosample_accessions'])
save('core_tree_metadata_complete.json',full);tsv('core_tree_metadata.tsv',wide);tsv('biosample_attributes_long.tsv',long)
save('publication_associations.json',associations);tsv('publication_associations.tsv',associations);save('doi_only_citations.json',doi_only)
used_pmids={a['pmid'] for a in associations}
curated=[dict(publications[p],associated_tree_tips=sorted({a['sample_id'] for a in associations if a['pmid']==p}),roles=sorted({a['role'] for a in associations if a['pmid']==p})) for p in sorted(used_pmids)]
save('curated_publications.json',curated);save('all_retrieved_publication_metadata.json',list(publications.values()))
tsv('publications.tsv',[{k:joined(v) if isinstance(v,list) else v for k,v in r.items() if k!='abstract'} for r in curated])
ris=[]
for c in curated:
    ris+=['TY  - JOUR','TI  - '+c['title'],'PY  - '+c['year'],'JO  - '+c['journal'],'DO  - '+c['doi'],'AN  - PMID:'+c['pmid'],'UR  - '+c['url']]
    ris+=['AU  - '+a for a in c['authors']];ris+=['N1  - Tree samples: '+joined(c['associated_tree_tips']),'N1  - Association roles: '+joined(c['roles']),'ER  - ','']
for c in doi_only:ris+=['TY  - JOUR','TI  - '+c['title'],'PY  - '+c['year'],'DO  - '+c['doi'],'UR  - '+c['url'],'N1  - '+c['evidence'],'ER  - ','']
(OUT/'curated_publications.ris').write_text('\n'.join(ris),encoding='utf-8')

group=collections.defaultdict(list)
for r in wide:
    if r['biosample_accessions']:group[r['biosample_accessions']].append(r['sample_id'])
duplicate_groups={k:v for k,v in group.items() if len(v)>1};save('shared_biosample_groups.json',duplicate_groups)
sentinels={'','missing','not collected','not provided','not applicable','unknown','na','n/a','none','not available','unspecified'}
def present(v):return bool(v) and any(x.strip().lower() not in sentinels for x in v.split('; '))
coverage={f:sum(present(r[f]) for r in wide) for f in fields}
summary={'tree_tips':66,'public_accession_mapped_tips':51,'local_unresolved_tips':15,'used_runs':29,'reference_assemblies':22,
         'unique_biosamples':len(group),'direct_bioprojects':len({p for r in wide for p in r['bioproject_accessions'].split('; ') if p}),
         'bioproject_records_downloaded':len(bp),'related_assembly_records':len(related['uids']),
         'related_sra_experiments':len({u for l in load('related_sra_links.json') for u in l['sra_uids']}),
         'curated_pubmed_publications':len(curated),'additional_doi_only_citations':len(doi_only),
         'tips_with_curated_publication_or_strain_background':len({a['sample_id'] for a in associations}),
         'shared_biosample_groups':duplicate_groups,'metadata_coverage_of_66_tips':coverage,
         'downloaded_files':len(load('download_manifest.json')),'downloaded_bytes':sum(d['bytes'] for d in load('download_manifest.json')),
         'scope':'All returned fields in linked NCBI records are retained. Search coverage is not a guarantee that every relevant publication has been indexed or found. Local-only identifiers were not queried externally.'}
save('summary.json',summary)

esc=html.escape
def link(url,label):return '<a href="'+esc(url,quote=True)+'">'+esc(label)+'</a>'
def accessions(s,kind):
    return '<br>'.join(link('https://www.ncbi.nlm.nih.gov/'+kind+'/'+a,a) for a in s.split('; ') if a)
def tbl(headers,rr,cls=''):
    return '<div class="scroll"><table class="'+cls+'"><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rr)+'</tbody></table></div>'
body='<header><p class="eyebrow">NCBI METADATA RETRIEVAL · 10 NOVEMBER 2026</p><h1>Core-tree sources and publications</h1><p>66 tree tips · 51 publicly accession-linked tips · 47 distinct BioSamples · 15 local isolates unresolved in NCBI</p></header>'
body+='<section><h2>What was retrieved</h2><p>Metadata were downloaded for the exact 66-tip core-tree cohort identified in the preceding audit. The 29 SRA run identifiers come directly from tree labels; 22 reference assembly accessions come from filenames in matching local source folders. This establishes an accession association, not a new proof that every historical annotation used exactly that genome file.</p>'
body+='<p>All returned fields are retained in raw XML/JSON and in the complete joined JSON. The catalogue includes Assembly/Datasets metadata, assembly reports/statistics, SRA experiment/study/sample/run metadata, BioSample attributes, BioProject records, taxonomic records, related assemblies/runs and publication metadata. Related public records are clearly separated from the inputs actually represented by tree tips.</p>'
body+=tbl(['Download','Contents'],[
 [link('core_tree_metadata.tsv','Metadata table (TSV)'),'66 rows; accessions, source/host/date/geography, sequencing fields, type-material statements, publication links and review flags.'],
 [link('core_tree_metadata_complete.json','Complete joined metadata (JSON)'),'All returned record fields joined to each tip, with local-source evidence and publication association roles.'],
 [link('biosample_attributes_long.tsv','All BioSample attributes (TSV)'),'Every returned attribute, its original name, sample association and NCBI source URL.'],
 [link('publications.tsv','Curated bibliography (TSV)')+' · '+link('curated_publications.ris','RIS for a reference manager'),f'{len(curated)} curated PubMed publications, plus one DOI-only historical strain paper in RIS/JSON.'],
 [link('publication_associations.tsv','Tip-to-publication mapping'),'Distinguishes genome publications, project associations, strain history, comparative reuse and taxonomic reassessment.'],
 [link('download_manifest.json','Download manifest'),'Source URLs, retrieval timestamps, byte counts and SHA-256 checksums for every downloaded response.'],
 [link('raw/','Raw download folder'),'Unmodified source responses, including unsuccessful/empty publication searches and contextual records.']])
body+='<p>No source project file was modified and no genome analysis was rerun. Metadata retrieval does not correct the unresolved alignment/annotation provenance gate from the earlier review.</p></section>'
body+='<section><h2>Findings that affect interpretation</h2><div class="alert"><b>Taxonomic conflict: IFO6365 and TN-484.</b> The 2024 study '+link('https://pubmed.ncbi.nlm.nih.gov/38445873/','Phylogenomics reveals extensive misidentification of fungal strains from the genus Aspergillus')+' explicitly places GCA_009932835.1 and GCA_009014675.2 under <i>A. pseudoterreus</i>, while the downloaded NCBI records still report <i>A. terreus</i>. Preserve both reported and literature-assessed labels. Validate the paired run-derived assemblies before changing all four related tree tips.</div>'
body+='<p><b>ASM-1 and ML-44 are biologically related.</b> The '+link('https://pubmed.ncbi.nlm.nih.gov/34485249/','ASM-1 study')+' identifies ASM-1 as a mutagenized derivative of ML-44. They should not automatically count as independent environmental samples in population-level comparisons.</p>'
body+='<p><b>Four assembly/run pairs share an NCBI BioSample.</b> This strengthens the identity mapping but does not make their assemblies or annotation versions interchangeable.</p>'
body+=tbl(['BioSample','Tree tips'],[[link('https://www.ncbi.nlm.nih.gov/biosample/'+k,k),esc(joined(v))] for k,v in duplicate_groups.items()])
body+='<p><b>Type-material fields are submitted metadata.</b> Nine mapped tips have type-material statements, including several FungID entries. These require strain/culture equivalence checks before using an ex-type symbol on a taxonomic tree. A pathogen-package label alone also does not prove that the isolate is clinically derived.</p>'
body+='<p><b>Publication linkage is graded.</b> Broad BioProject PRJNA176382 returns unrelated fungal-genome papers; these are retained in raw downloads but excluded from the curated bibliography. Geotrichum-only references to PRJNA833221 and a section Nigri paper mentioning PRJNA1265142 are not assigned to your terreus tips merely because they share a project. The gut-fungal project article is labelled project-associated until strain-specific supplementary confirmation.</p>'
body+='<p><b>The 15 local isolates remain unmapped:</b> '+esc(joined([r['sample_id'] for r in wide if not r['biosample_accessions']]))+'. No local-only DTO identifiers were sent to external services. This means public linkage is unresolved, not that the isolates are necessarily unpublished.</p></section>'
body+='<section><h2>Metadata available for figures</h2><p>Counts below are tree tips with an informative submitted value; alternative assemblies of the same BioSample count twice here. Placeholder values remain intact in downloads but are excluded from these coverage counts. Do not substitute NCBI submission dates for collection dates.</p>'
body+=tbl(['Field','Tips with informative values / 66'],[[esc(f),str(coverage[f])] for f in ['geo_loc_name','collection_date','isolation_source','host','host_disease','strain','isolate','type_material','ploidy','sample_type']])
body+='<p>Useful tree tracks are study/reference status, submitted isolation source, host category, country, date precision, sequencing technology and publication provenance. Keep raw submitted taxonomy distinct from reviewed taxonomy and type status. Clinical and environmental phenotypes from publications should be imported only after exact strain matching and assay comparability checks.</p></section>'
body+='<section><h2>All 66 tree tips</h2><p>Search by tip, accession, country, source, host or PMID. An empty field means no value was retrieved for that mapped record; original NCBI missing-value strings are preserved.</p><input id="filter" placeholder="Filter the metadata table…" aria-label="Filter metadata table">'
body+=tbl(['Tree tip','Accessions','Submitted identity / type','Geography / collection','Source / host','Publications / flags'],[
 [esc(r['sample_id']),accessions(r['used_run_accessions'],'sra')+accessions(r['local_reference_assembly_accessions'],'datasets/genome')+accessions(r['biosample_accessions'],'biosample')+accessions(r['bioproject_accessions'],'bioproject'),
  esc(r['ncbi_organism'])+'<br>'+esc(r['strain'] or r['isolate'])+'<br><small>'+esc(r['type_material'])+'</small>',
  esc(r['geo_loc_name'])+'<br>'+esc(r['collection_date']),esc(r['isolation_source'])+'<br>'+esc(r['host'])+'<br><small>'+esc(r['host_disease'])+'</small>',
  '<br>'.join(link('https://pubmed.ncbi.nlm.nih.gov/'+p+'/',p) for p in r['publication_pmids'].split('; ') if p)+'<br><small>'+esc(r['taxonomy_review_flag'] or r['biological_relationship_flag'] or ('Local public linkage unresolved' if not r['biosample_accessions'] else r['publication_roles']))+'</small>'] for r in wide],cls='sampletable')+'</section>'
body+='<section><h2>Curated publications</h2><p>Genome/phenotype publications support the named data or strains. Project-level associations are broader; background papers and comparative reuse must not be mistaken for the original genome source. The association table contains the evidence and caveats for each link.</p>'
body+=tbl(['Publication','Associated tree tips','Role'],[[link(c['url'],c['title'])+'<br><small>'+esc(c['year']+' · '+c['journal'])+'<br>'+link('https://doi.org/'+c['doi'],c['doi'])+'</small>',esc(joined(c['associated_tree_tips'])),esc(joined(c['roles']))] for c in sorted(curated,key=lambda x:(x['year'],x['title']),reverse=True)])
body+='<p>Additional historical citation: '+link(doi_only[0]['url'],doi_only[0]['title'])+' (1995). This describes TN-484 strain development, not genome sequencing.</p></section>'
body+='<section><h2>Remaining ambiguity and retrieval scope</h2><ul><li>Public metadata are linked for all 51 tips carrying established public sources. The 15 local-only tips require supplied accessions or separately authorized identifier searches.</li><li>Not every public tip has a located publication. A missing link means no verified association was found by this retrieval, not proof that no paper exists.</li><li>Some references have alternative assembly versions, shared BioSamples, or source/host/type annotations needing review. Exact historical annotation-to-genome hashes remain unresolved.</li><li>NCBI metadata, assembly statistics and taxonomy are a current snapshot. Archive release dates and annotation versions are retained; organelle inclusion and paired RefSeq/GenBank definitions can affect reported statistics.</li><li>Publication search used NCBI cross-links, exact public-accession searches in PubMed and PMC, and targeted title/strain follow-up. Raw candidate/context hits are retained; unrelated umbrella-project and cited-reference links were excluded from the curated set.</li><li>The downloaded package retains all fields returned by the queried record types. It does not claim exhaustive retrieval of every external laboratory record, supplementary phenotype table, sequence feature or publication ever associated with these strains.</li></ul></section>'
body+='<footer>Read-only source review. '+str(summary['downloaded_files'])+' downloaded responses, '+f"{summary['downloaded_bytes']/1024/1024:.1f}"+' MiB. Metadata and publication provenance are preserved in the download manifest.</footer>'
css='''body{font:15px/1.6 "Segoe UI",Arial,sans-serif;background:#f2f6f8;color:#203443;margin:0}main{max-width:1450px;margin:auto;padding:28px}header{background:#12384a;color:white;border-radius:15px;padding:32px}h1{font-size:38px;margin:12px 0}.eyebrow{letter-spacing:2px;font-size:12px;color:#90d7da}section{background:white;border:1px solid #dce5eb;border-radius:12px;padding:26px;margin:20px 0}h2{color:#174a5d}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:12px;vertical-align:top;border-bottom:1px solid #dbe5ea;text-align:left}th{background:#174a5d;color:white}tbody tr:nth-child(even){background:#f3f7f9}a{color:#00758b;overflow-wrap:anywhere}small{color:#536c7a}input{width:95%;padding:12px;border:1px solid #b6c9d2;border-radius:7px;margin:12px 0 22px;font:inherit}.alert{padding:20px;border-left:5px solid #c97828;background:#fff1dc;border-radius:6px}footer{color:#586d79;padding:20px}.sampletable td:first-child{font-weight:600;min-width:120px}.sampletable td{min-width:140px}@media print{input{display:none}body{background:white}section{border:0}table{font-size:9px}tr{break-inside:avoid}}'''
script="document.getElementById('filter').addEventListener('input',function(){const q=this.value.toLowerCase();document.querySelectorAll('.sampletable tbody tr').forEach(r=>r.style.display=r.textContent.toLowerCase().includes(q)?'':'none')});"
(OUT/'ncbi_metadata_report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Terreus core tree — NCBI metadata</title><style>'+css+'</style><main>'+body+'</main><script>'+script+'</script></html>',encoding='utf-8')

# Verify every raw response remains byte-identical to the retrieval manifest.
for d in load('download_manifest.json'):
    p=OUT/d['file'];assert p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==d['sha256'],p
save('validation.json',{'tree_tip_count':len(wide),'unique_tip_count':len({r['sample_id'] for r in wide}),
 'mapped_public_tips':sum(bool(r['biosample_accessions']) for r in wide),'raw_download_checksums':'all passed',
 'publication_associations_resolve':all(a['pmid'] in publications for a in associations),'source_project_modified':False})
archive=OUT.parent/'ncbi_core_tree_metadata_bundle.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.rglob('*'):
        if p.is_file():z.write(p,Path(OUT.name)/p.relative_to(OUT))
print(json.dumps(summary,indent=2,ensure_ascii=False))
print('REPORT',OUT/'ncbi_metadata_report.html')
print('ARCHIVE',archive)
