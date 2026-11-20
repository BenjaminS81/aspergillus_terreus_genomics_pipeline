from harvest_ncbi_core_metadata import *
hits=json.loads((OUT/'pmc_accession_search_hits.json').read_text(encoding='utf-8'))
uids=sorted({x for h in hits for x in h['pmc_uids']})
# Public full-text records are retained solely as source evidence for accession mentions.
root=eu('efetch',{'db':'pmc','id':','.join(uids),'retmode':'xml'},'pmc_accession_evidence.xml','xml')
evidence=[]
for a in root.findall('.//article'):
    pmcid=next((x.text for x in a.findall('./front/article-meta/article-id') if x.get('pub-id-type') in ['pmc','pmcid']),'')
    uid=(pmcid or '').removeprefix('PMC')
    relevant=[h for h in hits if uid in h['pmc_uids']]
    terms=sorted({t for h in relevant for t in h['public_accession_terms']})
    matches=[]
    for path,kind in [('./body//p','body'),('./body//table-wrap','table'),('./back//sec','back_section'),('./back//ref','reference'),('./back//notes','back_notes'),('./floats-group//table-wrap','float_table')]:
        for node in a.findall(path):
            s=' '.join(''.join(node.itertext()).split())
            matched=[t for t in terms if t in s]
            if matched:matches.append({'section_kind':kind,'accessions':matched,'text':s})
    evidence.append({'pmcid':'PMC'+uid,'title':text(a,'front/article-meta/title-group/article-title'),
                     'query_matches':relevant,'exact_accession_passages':matches,
                     'sample_mentions':{s:bool(re.search(r'(?<![\w])'+re.escape(s)+r'(?![\w])',' '.join(a.itertext()),re.I)) for h in relevant for s in h['samples']}})
save('publication_accession_evidence.json',evidence)
print('PUBLICATION EVIDENCE',len(evidence),flush=True)
for e in evidence:
    print(e['pmcid'],e['title'],[(m['section_kind'],m['accessions']) for m in e['exact_accession_passages']],flush=True)
    if e['pmcid'] in ['PMC13640352','PMC13610854','PMC9814666']:
        for m in e['exact_accession_passages']:print(m['text'][:1400],flush=True)
