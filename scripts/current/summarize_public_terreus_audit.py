"""Reconcile saved public metadata with read-only local evidence; no genome mutation."""
from pathlib import Path
import json,re,collections,html,datetime
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'review_reports/2026-11-20_public_genome_audit_v1'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
n=read(OUT/'ncbi_assemblies.json'); ena=read(OUT/'ena_assembly.json');runs=read(OUT/'ena_read_run.json')
files=read(OUT/'local_assembly_files.json');raw=read(OUT/'local_raw_read_files.json')
sources=read(ROOT/'taxonomy_reference_build/2026-11-10/analysis/genome_sources.json')
elig=read(ROOT/'gate_preparation/2026-11-18_v1/data/analysis_eligibility.json')
meta=read(ROOT/'review_reports/2026-11-10/ncbi_core_tree_metadata/core_tree_metadata_complete.json')
annotations=read(ROOT/'stakeholder_release/2026-11-18_v1_complete_three_gene/data/new_tree_annotations.json')
def runids(s):return set(re.findall(r'[SED]RR\d+',s))
def accs(s):return set(re.findall(r'GC[AF]_\d+\.\d+',s))
runmap=collections.defaultdict(set);bsmap=collections.defaultdict(set)
for s in sources:
 for r in runids(s['sample']):runmap[r].add(s['sample'])
for m in meta:
 s=m['summary']
 for b in re.findall(r'SAM[NED][A-Z]*\d+',s.get('biosample_accessions','')):bsmap[b].add(s['sample_id'])
for r in runs:bsmap[r['sample_accession']].update(runmap[r['run_accession']])
notes={
 'GCA_002930435.1':'Downloaded in inbox/strains/ncbi_dataset(4); absent from selected genome source manifest. Candidate for QC and inclusion.',
 'GCA_988135485.1':'Public release not found. BioSample maps to ERR17672364, an existing tree input. Its canonical sequence hash equals the ERR14883637 representative in the existing eligibility audit. Resolve assembly/read provenance before treating this as a new independent genome.',
 'GCA_061584345.1':'CMI 44339: no selected source linked by accession/BioSample; candidate new reference, subject to QC and taxonomic verification.',
 'GCA_986281685.1':'DAH1047FM: no selected source linked by accession/BioSample; candidate new reference, subject to QC and taxonomic verification.',
 'GCA_000166915.1':'Historical ATCC 20542 release. ATCC 20542 is already represented by GCA_016808415.1 and the SRR9691935-derived assembly. This is not an additional independent culture merely because its accession differs.'}
assemblies=[]
for x in n:
 a=x['accession']
 if a.startswith('GCF_'):continue
 aliases={a,x.get('paired_accession','')}-{''}
 f=[z for z in files if z['accession'] in aliases and z['nonempty']]
 selected=[s for s in sources if accs(s['windows_path'])&aliases]
 info=x['assembly_info'];bio=info.get('biosample',{}).get('accession','')
 attr={z['name']:z.get('value','') for z in info.get('biosample',{}).get('attributes',[])}
 infra=x['organism'].get('infraspecific_names',{})
 label=infra.get('strain') or infra.get('isolate') or attr.get('sample_name') or info.get('assembly_name','')
 assemblies.append({'accession':a,'paired_refseq':x.get('paired_accession',''),'label':label,'biosample':bio,'bioproject':info.get('bioproject_accession',''),'release_date':info.get('release_date',''),'assembly_level':info.get('assembly_level',''),'assembly_status':info.get('assembly_status',''),'available_exact_release':bool(f),'files':f,'selected_sources':[s['sample'] for s in selected],'in_66_tip_core':[s['sample'] for s in selected if s['in_core_tree']],'in_230_tip_three_gene':[s['sample'] for s in selected if 'LOCAL_'+s['sample'] in annotations],'biosample_linked_sources':sorted(bsmap[bio]),'status':'available_selected' if selected and f else 'available_not_selected' if f else 'release_not_found','notes':notes.get(a,''),'ncbi_url':'https://www.ncbi.nlm.nih.gov/datasets/genome/'+a+'/','ena_url':'https://www.ebi.ac.uk/ena/browser/view/'+a})
runrows=[]
for r in runs:
 rid=r['run_accession'];f=[z for z in raw if z['run_accession']==rid and z['bytes']>0]
 names={Path(z['path']).name for z in f};expected=[u.rsplit('/',1)[-1] for u in r.get('fastq_ftp','').split(';') if u]
 exact=bool(expected) and all(x in names for x in expected)
 samples=sorted(runmap[rid]);asm=[a['accession'] for a in assemblies if a['biosample']==r['sample_accession'] and a['available_exact_release']]
 note=''
 if r.get('instrument_platform')=='CAPILLARY':note='Legacy capillary run; assess study/sample grouping and aggregate coverage rather than counting every run as a separate modern genome.'
 if rid in ['SRR30681804','SRR30681805']:note='No.211 sector/wildtype Nanopore data: candidate new genomes after study design and assembly QC review.'
 if rid=='SRR9691936':note='Both paired FASTQs are present. ENA additionally lists a singleton FASTQ without _1/_2 suffix; that file was not found locally. This is not a missing paired mate.'
 if samples and not f:note+=' Selected assembly exists under this run label, but no accession-named raw FASTQ found; source/alias verification required.'
 runrows.append({**r,'local_files':f,'ena_fastq_filenames_all_present':exact,'selected_sources':samples,'same_biosample_available_assemblies':asm,'status':'raw_files_found' if f else 'assembly_source_present_raw_not_found' if samples else 'same_biosample_assembly_present_raw_not_found' if asm else 'raw_not_found_no_selected_source_link','notes':note.strip()})
atcc=[]
for tip in ['LOCAL_ATCC20541','LOCAL_ATCC20542','LOCAL_SRR9691935','LOCAL_SRR25621424','LOCAL_SRR25621436','LOCAL_SRR25636132','LOCAL_SRR25636149']:
 ann=annotations[tip];sample=tip.removeprefix('LOCAL_');r=next((r for r in runs if r['run_accession']==sample),{})
 atcc.append({'tree_tip':tip,'reported_culture':r.get('strain') or sample,'in_core_source_manifest':next(s['in_core_tree'] for s in sources if s['sample']==sample),'in_complete_three_gene_tree':True,'three_gene_type_clade':ann['new_clade_category'],'clade_code':ann['new_clade_code'],'formal_species_assignment':ann['species_ring'],'provenance':'Public assembly' if sample.startswith('ATCC') else 'Local assembly from public reads; not proof of laboratory-generated sequencing'})
ns={x['accession'] for x in n};es={x['assembly_accession']+'.'+x['version'] for x in ena}
summary={'audit_date':'2026-11-20','scope':'Public records assigned to NCBI taxon 33178 (Aspergillus terreus) or its descendants, plus ENA WGS runs under that taxonomy subtree. Does not assert completeness for all section Terrei species or reclassified records outside this subtree.','ncbi_dataset_records':len(n),'ncbi_distinct_genbank_assemblies':len(assemblies),'ena_assembly_records':len(ena),'ncbi_only_accessions':sorted(ns-es),'ena_only_accessions':sorted(es-ns),'available_current_assembly_accessions':sum(a['available_exact_release'] for a in assemblies),'selected_public_assembly_accessions':sum(bool(a['selected_sources']) for a in assemblies),'available_not_selected':[a['accession'] for a in assemblies if a['status']=='available_not_selected'],'public_releases_not_found':[a['accession'] for a in assemblies if not a['available_exact_release']],'ena_wgs_runs':len(runs),'ena_wgs_biosamples':len({r['sample_accession'] for r in runs}),'ena_wgs_run_status_counts':dict(collections.Counter(r['status'] for r in runrows)),'ena_runs_with_all_expected_fastq_names':sum(r['ena_fastq_filenames_all_present'] for r in runrows),'legacy_capillary_runs':sum(r.get('instrument_platform')=='CAPILLARY' for r in runs),'local_scan':read(OUT/'local_scan.json'),'limits':['Nonempty files and accession/filename evidence establish availability, not verified sequence provenance or file integrity. No new genome/read downloads or computational jobs launched.','Selected source manifest contains 77 labels across the broader study; exact sequence deduplication, independent cultures and tree-specific eligibility are separate counts.','NCBI includes the paired NIH2624 RefSeq record; ENA instead also includes superseded TN-484 version 1. Neither is an extra independent current assembly.','ATCC is a culture collection identifier, not a species determination or evidence that sequencing was performed in the user laboratory.','ATCC 1012 is reported for four FungID runs that occupy different marker clades; retain as unresolved source-identity/taxonomy evidence, not proven interchangeable cultures.']}
save('assembly_crosscheck.json',assemblies);save('wgs_run_crosscheck.json',runrows);save('atcc_tree_crosscheck.json',atcc);save('audit_summary.json',summary)
def esc(x):return html.escape(str(x))
def table(headers,rows):return '<table><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in r)+'</tr>' for r in rows)+'</tbody></table>'
page='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Terreus public genome audit — 2026-11-20</title><style>body{font:15px/1.5 system-ui,sans-serif;max-width:1450px;margin:35px auto;padding:0 25px;color:#172d40}h1,h2{color:#164d62}table{border-collapse:collapse;width:100%;font-size:13px;margin:18px 0}th,td{padding:9px;border:1px solid #c9d4dd;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#e7f0f4}tr:nth-child(even){background:#f6f8fa}.lead{background:#eaf2f5;padding:18px}a{color:#006890}small{color:#536571}</style>'''
page+='<h1>Aspergillus terreus: public genome availability audit</h1><p>Reviewed 20 November 2026. Original project data and active scientific jobs were not modified.</p>'
page+='<p class="lead"><b>24 of 28 distinct current GenBank assembly accessions were found locally.</b> Of these, 23 are selected public assembly sources. T3_Kankrej is downloaded in the inbox but not selected. Four published releases were not found; they are not necessarily four additional independent strains.</p>'
page+='<h2>ATCC identity and taxonomic placement</h2><p>ATCC 20541 and ATCC 20542 are public records of ATCC-labelled cultures and are present in both tree datasets. The current three-gene placements are respectively the A. montoensis and A. citrinoterreus type-containing clades; formal species assignment remains unresolved. SRR9691935 is an additional public-read representation of ATCC 20542. No sample key establishes an independently sequenced laboratory counterpart.</p>'
page+=table(['Tree entry','Reported culture','Core source','Complete three-gene tree','Current three-gene type clade','Provenance'],[(a['tree_tip'],a['reported_culture'],a['in_core_source_manifest'],a['in_complete_three_gene_tree'],a['three_gene_type_clade'],a['provenance']) for a in atcc])
page+='<p><b>Additional identity gate:</b> Four FungID entries report ATCC 1012 but occupy different marker clades. Shared culture text is insufficient to collapse these records or conclude that all four represent a verified identical culture.</p>'
page+='<h2>Assembly gaps and unselected downloads</h2>'+table(['Accession','Reported strain/sample','Availability','BioSample-linked selected source','Interpretation'],[(a['accession'],a['label'],a['status'],', '.join(a['biosample_linked_sources']),a['notes']) for a in assemblies if a['status']!='available_selected'])
page+='<h2>Complete public assembly cross-check</h2>'+table(['Accession','Reported strain/sample','BioSample','Release','Available','Selected source','Core source / complete 3-gene source'],[(a['accession'],a['label'],a['biosample'],a['release_date'],a['available_exact_release'],', '.join(a['selected_sources']),', '.join(a['in_66_tip_core'])+' / '+', '.join(a['in_230_tip_three_gene'])) for a in assemblies])
page+='<h2>NCBI / ENA reconciliation</h2><p>NCBI Datasets returned 29 records: 28 GenBank assemblies plus the paired NIH2624 RefSeq record GCF_000149615.1. ENA returned the same 28 GenBank accessions plus the superseded TN-484 version GCA_009014675.1. The locally selected TN-484 release is version 2. No additional unique current assembly was supplied by ENA in this taxonomy query.</p>'
page+='<h2>Public WGS reads</h2><p>'+esc(f"ENA returned {len(runs)} runs across {summary['ena_wgs_biosamples']} BioSamples. {sum(bool(r['local_files']) for r in runrows)} runs have accession-named raw FASTQ files in the scanned local folders; {sum(r['ena_fastq_filenames_all_present'] for r in runrows)} have all expected ENA filenames. Missing raw files do not imply a missing genome when an assembly is present. Of the returned runs, 36 are legacy capillary sequencing runs, not 36 independent new genomes.")+'</p>'
page+=table(['Run','BioSample','Reported strain','Platform','Local status','Selected source','Notes'],[(r['run_accession'],r['sample_accession'],r.get('strain',''),r.get('instrument_platform',''),r['status'],', '.join(r['selected_sources']),r['notes']) for r in runrows])
page+='<h2>Recommended next gates</h2><ol><li>Evaluate T3_Kankrej, already downloaded, for genome/locus QC and tree inclusion.</li><li>Obtain and QC CMI 44339 and DAH1047FM if expanding the cohort.</li><li>Resolve 59TF / ERR17672364 versus T59-F / ERR14883637 provenance before adding GCA_988135485.1 as an independent genome. The existing eligibility audit reports identical canonical sequence hashes; ENA compressed FASTQ checksums differ, so exact read equivalence is not established.</li><li>Treat historical ATCC 20542 GCA_000166915.1 as a reference-comparison option, not automatically an additional strain.</li><li>Review the No.211 wildtype/sector Nanopore pair (SRR30681805 / SRR30681804) as potential additional genomes, keeping the experimental relationship explicit.</li><li>Recover or document raw-read aliases for DRR191518, DRR204676 and ERR17672364 where read-based QC is needed; their assembly source labels already exist.</li><li>Preserve reported ATCC identities separately from proposed species placement and resolve the discordant ATCC 1012 metadata.</li></ol>'
page+='<h2>Scope and evidence</h2><p>'+esc(summary['scope'])+'</p><p>'+esc(summary['local_scan']['scope'])+'</p><ul>'+''.join('<li>'+esc(s)+'</li>' for s in summary['limits'])+'</ul>'
page+='<p>Machine-readable evidence: <a href="audit_summary.json">summary</a> · <a href="assembly_crosscheck.json">assemblies and local paths</a> · <a href="wgs_run_crosscheck.json">all runs</a> · <a href="atcc_tree_crosscheck.json">ATCC tree entries</a> · <a href="retrieval_manifest.json">API queries, timestamps and SHA256 checksums</a>.</p>'
page+='<p>Public sources: <a href="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/33178/dataset_report?page_size=1000">NCBI Datasets taxon 33178</a> · <a href="https://www.ebi.ac.uk/ena/portal/api/">ENA Portal API</a> · <a href="https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/taxon-based-search.html">ENA taxon search documentation</a>. Exact queries and raw responses are retained in this folder.</p></html>'
(OUT/'public_genome_audit.html').write_text(page,encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k not in ['limits','local_scan','scope']},indent=2))
