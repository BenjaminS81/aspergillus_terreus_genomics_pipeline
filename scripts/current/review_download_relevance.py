"""Review downloaded data against actual released tree tips and public QC metadata."""
from pathlib import Path
import json,re,hashlib,html,datetime,collections
R=Path(__file__).resolve().parent
OLD=R/'review_reports/2026-11-20_public_genome_audit_v1'
OUT=R/'review_reports/2026-11-20_public_genome_audit_v2';OUT.mkdir(exist_ok=True)
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
treepaths={'core66':R/'stakeholder_release/2026-11-17_v1/data/core_genome_original.treefile','complete230':R/'stakeholder_release/2026-11-18_v1_complete_three_gene/runs/complete230.treefile','full283':R/'stakeholder_release/2026-11-18_v1_complete_three_gene/runs/full283_existing/combined_all_strains.treefile'}
tips={k:set(re.findall(r'(?:\(|,)([^():,;]+):',p.read_text())) for k,p in treepaths.items()}
assert {k:len(v) for k,v in tips.items()}=={'core66':66,'complete230':230,'full283':283}
tree_evidence=[{'tree':k,'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'tip_count':len(tips[k])} for k,p in treepaths.items()]
sources=load(R/'taxonomy_reference_build/2026-11-10/analysis/genome_sources.json');n=load(OLD/'ncbi_assemblies.json')
assemblies=load(OLD/'assembly_crosscheck.json');raw=load(OLD/'local_raw_read_files.json');runs=load(OLD/'wgs_run_crosscheck.json')
extra=load(OUT/'additional_downloaded_runs_ena.json');ann=load(R/'stakeholder_release/2026-11-18_v1_complete_three_gene/data/new_tree_annotations.json')
rows=[]
for s in sources:
 name=s['sample'];pub=s['source_kind']=='public_reference_genome' or bool(re.search(r'[SED]RR\d+',name));f=Path(s['windows_path'])
 rows.append({'sample':name,'public_data':pub,'provenance':s['source_kind'],'run_accessions':re.findall(r'[SED]RR\d+',name),'assembly_accessions':re.findall(r'GC[AF]_\d+\.\d+',s['windows_path']),'assembly_path':str(f),'assembly_present':f.is_file() and f.stat().st_size>0,'core66':name+'_transcripts' in tips['core66'],'complete230':'LOCAL_'+name in tips['complete230'],'full283':'LOCAL_'+name in tips['full283']})
recommendations={
 'GCA_061584345.1':('High: download assembly','CMI 44339 is a genome-sized, reasonably contiguous assembly (29.88 Mb; 153 scaffolds; scaffold N50 0.886 Mb). It adds a named culture and a connection to experimental terramide biosynthesis. Run standard genome/locus QC before inclusion; the terramide study does not itself establish a lovastatin phenotype.','https://www.nature.com/articles/s42004-024-01311-2.pdf'),
 'GCA_988135485.1':('High: download for provenance reconciliation','59TF is genome-sized (29.99 Mb, 179 contigs). Its BioSample links to ERR17672364 already in the trees. The locally assigned ERR17672364 assembly shares a canonical sequence hash with ERR14883637; compare the published releases before treating this as an independent genome. The two public assemblies have different lengths/contig statistics.',''),
 'GCA_986281685.1':('Do not add as a genome','DAH1047FM is only 347 bp, sequenced by ABI 3500xl and assembled with Chromas; NCBI flags genome length too small. Despite the submitted Complete Genome label, this is not a usable whole-genome reference. Locus identity would need separate sequence verification before any marker use.',''),
 'GCA_000166915.1':('Low: historical comparison only','Older ATCC 20542 release: only 138,520 bp across 175 contigs, NCBI partial flag. The culture is already represented by a modern public assembly and public-read assembly. It adds neither a complete genome nor an independent strain.',''),
 'GCA_002930435.1':('Already downloaded: conditional QC / marker recovery','T3_Kankrej is in inbox but absent from both released tree datasets. Highly fragmented: 13,340 contigs, contig N50 3,497 bp, reported coverage 6.87x. Attempt validated marker recovery and assess completeness/contamination before any genome-tree inclusion. No re-download needed.',''),
 'GCA_056506665.1':('Already downloaded: retain exclusion pending marker rescue','HNGD-TM15 is present in the broader 283-tip tree but absent from core66 and complete230. The selection audit contains BenA only (CaM and RPB2 absent). Inclusion in a complete three-locus tree requires validated recovery, not simply a new download.','')}
for a in assemblies:
 record=next(x for x in n if x['accession']==a['accession']);a['assembly_stats']=record['assembly_stats'];a['ncbi_warnings']=record['assembly_info'].get('genome_notes',[])
 linked=[r for r in rows if a['accession'] in r['assembly_accessions']]
 a['actual_core66_sources']=[r['sample'] for r in linked if r['core66']];a['actual_complete230_sources']=[r['sample'] for r in linked if r['complete230']];a['actual_full283_sources']=[r['sample'] for r in linked if r['full283']]
 a['recommendation'],a['rationale'],a['publication']=recommendations.get(a['accession'],('Already included: no assembly download needed','Retain current selected assembly, subject to the ongoing QC and provenance gates.',''))
 # Replace preliminary v1 notes where the quality review supersedes them.
 if a['accession'] in recommendations:a['notes']=a['rationale']
runrows=[]
for r in runs+extra:
 rid=r['run_accession'];linked=[s for s in rows if rid in s['run_accessions']];local=[x for x in raw if x['run_accession']==rid]
 x={**r,'local_files':local,'actual_core66_sources':[s['sample'] for s in linked if s['core66']],'actual_complete230_sources':[s['sample'] for s in linked if s['complete230']],'actual_full283_sources':[s['sample'] for s in linked if s['full283']]}
 x['download_recommendation']='Already downloaded; no re-download indicated' if local else 'Assess sample/study relationship before requesting additional data'
 if rid in ['DRR191518','DRR204676']:x['download_recommendation']='Optional: recover reads for read-backed QC; corresponding public assemblies and local run-labelled assemblies already included. Does not add a new strain.'
 if rid=='ERR17672364':x['download_recommendation']='Conditional: prioritize assembly provenance reconciliation; recover reads if needed to distinguish 59TF from T59-F. About 0.80 GB compressed.'
 if rid in ['SRR30681804','SRR30681805']:x['download_recommendation']='High scientific relevance, conditional download: Terrnet 211 wildtype/sector pair from Eisele et al. 2025; useful phenotype/secondary-metabolism contrast. Check original study archives before re-downloading. Approximately 7.00 GB compressed combined. Keep as a related pair, not two independent taxon samples.'
 if r.get('instrument_platform')=='CAPILLARY':x['download_recommendation']='Low priority: legacy capillary runs grouped across two BioSamples; not 36 additional genomes. NIH2624 is already represented. Resolve the other BioSample before assuming it adds a new strain.'
 if rid=='SRR9691936':x['download_recommendation']='No substantial data gap: both paired files are present. ENA additionally lists a 456-byte compressed singleton file; negligible for genome expansion, optional for exact archival completeness.'
 if r in extra:
  x['download_recommendation']='Already downloaded: optional broader comparison/outgroup assembly; no re-download. No run-derived genome appears in current trees.'
  if rid=='SRR10511972':x['taxon_caveat']='ENA run scientific_name/tax_id fields are blank; A. neoindicus CBS 444.75 is the study title, pending sample-level confirmation.'
  norm=lambda s:re.sub(r'[^a-z0-9]','',s.lower())
  culture=r.get('strain') or ('CBS 444.75' if rid=='SRR10511972' else '')
  x['literature_reference_tips_same_culture_text']=[k for k,v in ann.items() if culture and norm(v.get('strain',''))==norm(culture) and not k.startswith('LOCAL_')]
 runrows.append(x)
save('actual_tree_inclusion_all_sources.json',rows);save('assembly_relevance_crosscheck.json',assemblies);save('run_relevance_crosscheck.json',runrows);save('tree_file_evidence.json',tree_evidence)
# One streaming pass through the already downloaded T3 assembly; no alignment or heavy QC job.
t3=next(a for a in assemblies if a['accession']=='GCA_002930435.1');f=Path(t3['files'][0]['path']);lengths=[];current=None
with f.open('rb') as handle:
 for line in handle:
  if line.startswith(b'>'):
   if current is not None:lengths.append(current)
   current=0
  else:
   assert current is not None
   current+=len(line.strip())
if current is not None:lengths.append(current)
total=sum(lengths);cumulative=0
for length in sorted(lengths,reverse=True):
 cumulative+=length
 if cumulative>=total/2:n50=length;break
save('t3_local_fasta_check.json',{'path':str(f),'records':len(lengths),'bases':total,'N50':n50,'matches_ncbi_contig_count_and_total':len(lengths)==int(t3['assembly_stats']['number_of_contigs']) and total==int(t3['assembly_stats']['total_sequence_length'])})
summary={'review_date':'2026-11-20','distinct_public_assembly_accessions_in_taxon_query':len(assemblies),'downloaded_assembly_accessions':sum(a['available_exact_release'] for a in assemblies),'downloaded_assemblies_in_core66':sum(bool(a['actual_core66_sources']) for a in assemblies),'downloaded_assemblies_in_complete230':sum(bool(a['actual_complete230_sources']) for a in assemblies),'public_source_labels':sum(x['public_data'] for x in rows),'public_source_labels_in_core66':sum(x['public_data'] and x['core66'] for x in rows),'public_source_labels_in_complete230':sum(x['public_data'] and x['complete230'] for x in rows),'local_public_raw_run_ids':len({x['run_accession'] for x in raw}),'local_public_raw_runs_without_selected_genome':['SRR10511972','SRR35106274','SRR6447866','SRR8397341','SRR8399221'],'correction_to_v1':'DAH1047FM is a 347-bp atypical record and is unsuitable for whole-genome addition. The historical ATCC20542 release is only 138.5 kb and partial. T3_Kankrej needs stringent QC because of severe fragmentation.','tree_inputs':tree_evidence,'scope':'Current NCBI/ENA A. terreus taxon subtree plus all public-labelled sources and accession-named downloaded raw reads identified in the local project audit. Not a new exhaustive public inventory of all section Terrei species.','no_new_genome_or_read_downloads':True,'scientific_jobs_unchanged':True}
save('relevance_summary.json',summary)
def e(x):return html.escape(str(x))
def table(head,body):return '<table><thead><tr>'+''.join('<th>'+e(x)+'</th>' for x in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+e(x)+'</td>' for x in row)+'</tr>' for row in body)+'</tbody></table>'
page='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Terreus download and inclusion review v2</title><style>body{font:15px/1.5 system-ui,sans-serif;max-width:1450px;margin:35px auto;padding:0 25px;color:#193348}h1,h2{color:#14536b}table{border-collapse:collapse;width:100%;font-size:13px;margin:18px 0}th,td{padding:9px;border:1px solid #c9d4dd;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#e7f0f4}tr:nth-child(even){background:#f6f8fa}.lead{background:#eaf2f5;padding:18px}a{color:#006890}</style>'''
page+='<h1>Downloaded data, actual tree inclusion and download relevance</h1><p>Version 2 — 20 November 2026. Supersedes the preliminary relevance recommendations in version 1. Genome/read downloads and running scientific jobs were unchanged.</p>'
page+='<p class="lead"><b>24 of 28 distinct public assembly accessions are downloaded; 22 appear in both released core66 and complete230 trees.</b> HNGD-TM15 is downloaded but lacks a complete marker set. T3_Kankrej is downloaded in the inbox but unselected and severely fragmented. Of four missing assembly releases, prioritize CMI 44339 for expansion and 59TF for provenance reconciliation. Do not treat the 347-bp DAH1047FM record as a genome.</p>'
page+='<h2>Decisions for missing or excluded assembly releases</h2>'+table(['Accession','Culture/sample','Core66 / complete230 / full283','Recommendation','Evidence'],[(a['accession'],a['label'],str(bool(a['actual_core66_sources']))+' / '+str(bool(a['actual_complete230_sources']))+' / '+str(bool(a['actual_full283_sources'])),a['recommendation'],a['rationale']) for a in assemblies if a['accession'] in recommendations])
page+='<p>T3_Kankrej local FASTA verification: '+e(f'{len(lengths):,} records; {total:,} bases; N50 {n50:,} bp.')+' The publication connection for CMI 44339 is <a href="https://www.nature.com/articles/s42004-024-01311-2.pdf">Han et al. 2024, terramide biosynthesis</a>; this is useful secondary-metabolism context, not evidence of its lovastatin phenotype.</p>'
page+='<h2>Additional downloaded raw data not represented by genome-derived tree tips</h2>'+table(['Run','Reported taxon / study title','Strain','Literature tip already in complete230','Action'],[(r['run_accession'],r.get('scientific_name') or r['study_title'],r.get('strain') or 'study title: CBS 444.75',', '.join(r.get('literature_reference_tips_same_culture_text',[])) or 'No exact culture-text match',r['download_recommendation']) for r in runrows if r['run_accession'] in summary['local_public_raw_runs_without_selected_genome']])
page+='<p>Four of these cultures already have literature-derived marker reference tips; those references do not mean the downloaded genome reads were assembled or included. Exact culture matching here is based on text and is not a new authentication claim. For SRR10511972, sample-level taxon fields are blank and the taxon/culture above comes from the study title.</p>'
page+='<h2>Raw-read downloads worth considering</h2><p>The missing Terrnet 211 wildtype/sector Nanopore pair (SRR30681805 / SRR30681804; about 7.00 GB compressed) has direct relevance to phenotype and secondary-metabolism comparisons. The pair belongs to <a href="https://journals.asm.org/doi/10.1128/mbio.03926-24">Eisele et al. 2025</a>. Check original study archives first; absence from this project does not prove absence elsewhere. Keep their experimental relationship explicit. RNA-seq is a separate follow-up for expression questions, not required for the genome/three-locus trees.</p><p>DRR191518, DRR204676 and ERR17672364 are assembly-source labels in the current trees but have no accession-named raw FASTQs in the scanned folders. Recovering reads improves provenance/read-backed QC but adds no established new strain. Legacy capillary runs are low priority. SRR9691936 has both paired files; its missing ENA singleton is only 456 compressed bytes.</p>'
page+='<h2>All 28 public assembly accessions</h2>'+table(['Accession','Culture/sample','Downloaded','Core66','Complete230','Full283','Recommendation'],[(a['accession'],a['label'],a['available_exact_release'],bool(a['actual_core66_sources']),bool(a['actual_complete230_sources']),bool(a['actual_full283_sources']),a['recommendation']) for a in assemblies])
page+='<h2>All public-labelled selected genome sources</h2>'+table(['Source label','Provenance','Assembly present','Core66','Complete230','Full283'],[(s['sample'],s['provenance'],s['assembly_present'],s['core66'],s['complete230'],s['full283']) for s in rows if s['public_data']])
page+='<p>These are source labels, not counts of independent cultures. Aliases, public assembly/read representations and identical sequences must remain distinguished. The 230-tip tree includes 75 genome-derived tips plus 155 literature/reference tips; downloaded assembly inclusion is a separate count.</p>'
page+='<h2>Scope and reproducibility</h2><p>'+e(summary['scope'])+'</p><p>Inclusion was checked against leaf labels in the actual saved tree files, not only source manifests. Local availability uses today’s scoped project inventory plus current existence checks of selected assemblies. It does not prove FASTQ integrity, sequence provenance, or absence of renamed files outside the inventoried locations. A BUSCO job remains active, so these counts describe the released tree files explicitly listed below.</p>'
page+=table(['Tree','Tips','File'],[(x['tree'],x['tip_count'],x['path']) for x in tree_evidence])
page+='<p>Evidence: <a href="relevance_summary.json">summary</a> · <a href="assembly_relevance_crosscheck.json">all assemblies, parameters and paths</a> · <a href="actual_tree_inclusion_all_sources.json">all source/tree memberships</a> · <a href="run_relevance_crosscheck.json">all reviewed raw runs</a> · <a href="tree_file_evidence.json">tree hashes</a> · <a href="../2026-11-20_public_genome_audit_v1/retrieval_manifest.json">original public API retrievals</a> · <a href="additional_retrieval_manifest.json">additional ENA retrieval</a>.</p></html>'
(OUT/'download_relevance_review.html').write_text(page,encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k!='tree_inputs'},indent=2))
