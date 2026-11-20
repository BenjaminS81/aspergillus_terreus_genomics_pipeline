from pathlib import Path
import json,csv,collections,hashlib
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10';A=P/'analysis';L=P/'reference_library'
def js(path):return json.loads(path.read_text())
def tsv(path,rows):
    if not rows:return
    cols=list(dict.fromkeys(k for d in rows for k in d))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
calls=js(A/'genomic_markers/locus_calls.json');mapped={(d['sample'],d['locus']):d for d in js(A/'orthology_crosscheck/validated_marker_gene_mapping.json')};old=js(A/'orthology_crosscheck/old_marker_genes.json');actual={d['sample']:d for d in js(A/'orthology_crosscheck/actual_old_BenA_alignment_mapping.json')}
rows=[]
for d in calls:
    m=mapped.get((d['sample'],d['locus']),{});g=m.get('genes',[]);o=[x for x in old if x['sample']==d['sample'] and x['locus']==d['locus']]
    if d['locus']=='BenA':o=[x for x in o if x['original_header'] in actual.get(d['sample'],{}).get('exact_raw_headers',[])]
    rows.append({'sample':d['sample'],'locus':d['locus'],'in_core_tree':d['in_core_tree'],'genomic_call':d['status'],'strong_genomic_candidates':d['strong_candidate_count'],'contig':d.get('contig',''),'start_1based':d.get('start_1based',''),'end_1based':d.get('end_1based_inclusive',''),'strand':d.get('strand',''),'query_coverage_percent':d.get('query_coverage_percent',''),'identity_percent':d.get('weighted_identity_percent',''),'annotation_status':m.get('status','not_in_annotation_mapping'),'validated_gene_ids':','.join(x['gene_id'] for x in g),'validated_orthogroups':','.join(sorted({str(x['orthogroup']) for x in g})),'old_extracted_gene_ids':','.join(x['old_gene'] for x in o),'old_extracted_orthogroups':','.join(sorted({str(x['orthogroup']) for x in o})),'genome_source':d['genome_source'],'genome_sha256':d['source_sha256']})
tsv(A/'marker_gene_crosscheck.tsv',rows)
groups=js(A/'assembly_sequence_groups.json');dups=[{'sequence_sha256':h,'labels':','.join(d['samples']),'file_count':len(d['files']),'searched_as':d['searched_as']} for h,d in groups.items() if len(d['samples'])>1]
tsv(A/'identical_assembly_label_groups.tsv',dups)
tsv(A/'assembly_asset_inventory.tsv',[{k:v for k,v in d.items() if k!='tag'} for d in js(A/'assembly_variant_inventory.json')])
refs=js(A/'tree_reference_metadata.json');reference_table=[]
for d in refs:
    row={'id':d['reference_id'],'analysis_label':d['analysis_species'],'strain':d['strain'],'role':d['role'],'type_flag':d['type_flag'],'type_name_in_accepted_table':d.get('accepted_type_name',''),'historical_type_name':d.get('historical_type_name',''),'section':d['section'],'series':d['series'],'ICPA_source_row':d.get('source_row',''),'historical_taxonomy_status':d.get('taxonomy_status','')}
    for l in ['BenA','CaM','RPB2']:row[l+'_accession']=d.get(l+'GB','');row[l+'_bases']=len(d[l+'_analysis_sequence'])
    reference_table.append(row)
tsv(L/'tree_reference_metadata.tsv',reference_table)
gb=js(L/'genbank_records.json');pubs=collections.defaultdict(set)
for a,d in gb.items():
    for ref in d['references']:
        if ref.get('title') and ref['title'].lower()!='direct submission':pubs[(ref.get('title',''),ref.get('journal',''),ref.get('pubmed',''))].add(a)
tsv(L/'marker_associated_publications.tsv',[{'title':title,'journal':journal,'PMID':pmid,'url':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/' if pmid else '', 'marker_accessions':','.join(sorted(acc))} for (title,journal,pmid),acc in pubs.items()])
(L/'downloaded_GenBank_marker_sequences.fasta').write_text(''.join('>'+d['accession-version']+' '+d['definition']+'\n'+d['sequence'].upper()+'\n' for d in gb.values()))
resources={'compute_cgroup':{'MemoryMax_bytes':2147483648,'MemoryHigh_bytes':1073741824,'MemorySwapMax_bytes':0,'CPUQuota_percent':100,'Nice':15,'IOWeight':10},'per_process_virtual_address_limit_bytes':1610612736,'concurrent_compute_workers':1,'source_project_policy':'read-only','Panaroo_policy':'No signals, settings, limits or files changed','iqtree_partition_memory_note':'IQ-TREE 2.4.0 rejects -mem with partition models; OS cgroup and RLIMIT_AS limits enforce the budget for the combined analysis.'}
(P/'resource_policy.json').write_text(json.dumps(resources,indent=2))
print('EXPORTED',len(rows),'locus checks;',len(reference_table),'tree reference records;',len(pubs),'publication entries; shared assembly groups',len(dups))
