from pathlib import Path
import sys,json,hashlib,subprocess,ast,zipfile,collections
OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
def load(p):return json.loads(p.read_text())
def save(p,v):p.write_text(json.dumps(v,indent=2),encoding='utf-8')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
summary=load(OUT/'reports/preparation_validation.json')
manifest=load(OUT/'data/analysis_eligibility.json')
cache=load(OUT/'data/current_input_hashes.json')
checks={}
# Verify every label path, including the five non-representative aliases.
for m in manifest:
 p=Path(m['assembly']);key=str(p.resolve()).lower();facts=cache.get(key)
 if not facts:
  seqs=[];h=hashlib.sha256();n=0
  for line in p.read_text().splitlines():
   if line.startswith('>'):
    if n:seqs.append(h.hexdigest())
    h=hashlib.sha256();n+=1
   else:h.update(line.strip().upper().encode())
  if n:seqs.append(h.hexdigest())
  facts={'file_sha256':sha(p),'sequence_sha256':hashlib.sha256('\n'.join(sorted(seqs)).encode()).hexdigest()}
 assert facts['sequence_sha256']==m['sequence_sha256']
 m['assembly_file_sha256']=facts['file_sha256']
save(OUT/'data/analysis_eligibility.json',manifest)
checks['all_77_label_assembly_hashes_verified']=len(manifest)==77
checks['all_72_sequence_groups_verified']=summary['canonical_sequence_groups']==72
checks['all_66_core_tree_labels_present']=sum(m['in_core_tree'] for m in manifest)==66
checks['no_busco_parse_errors']=summary['busco_parse_errors']==0
checks['partition_complete']=summary['partition_blocks']==4661 and summary['partition_columns']==9012257
for p in (OUT/'scripts').glob('*.py'):ast.parse(p.read_text())
checks['python_syntax']=True
import tempfile
tests=Path(tempfile.mkdtemp(prefix='primary_selector_test_',dir=OUT/'reports'))
(tests/'input.faa').write_text('>a1\nMK\n>a2\nMKK\n>b1\nMKK\n>c1\nMK*K\n')
(tests/'input.gff3').write_text('##gff-version 3\n'+''.join(f'c\tt\tmRNA\t1\t9\t.\t+\t.\tID={t};Parent={g}\n' for t,g in [('a1','g1'),('a2','g1'),('b1','g2'),('c1','g3')]))
result=subprocess.run([sys.executable,str(OUT/'scripts/export_primary_proteins.py'),'--protein',str(tests/'input.faa'),'--gff',str(tests/'input.gff3'),'--output',str(tests/'out')],capture_output=True,text=True,check=True)
sel=load(tests/'out/selection.json')
checks['longest_isoform_selected']=any(r['protein_id']=='a2' and r['selected'] for r in sel['mapping'])
checks['identical_paralogs_not_collapsed']=sel['selected_genes']==2
checks['internal_stop_rejected']=len(sel['rejected'])==1 and sel['rejected'][0]['reason']=='internal_stop'
queue=load(OUT/'data/busco_missing_genome_queue.json')
if queue:
 r=subprocess.run([sys.executable,str(OUT/'scripts/run_one_busco.py'),queue[0]['job']],capture_output=True,text=True,check=True)
 checks['busco_default_dry_run']=not (OUT/'runs').exists() and queue[0]['job'] in r.stdout
checks['no_analyses_launched']=summary['analyses_launched']==0 and not (OUT/'runs').exists()
assert all(checks.values()),checks
save(OUT/'reports/package_validation.json',{'checks':checks,'resource_limits':LIMITS,'heavy_analysis_and_runtime_cgroup_test':'not executed'})

annotation=load(OUT/'data/annotation_candidates.json');proposals=[]
canonical={m['sample']:m['sequence_sha256'] for m in manifest}
def seqhash(path):
 key=str(path.resolve()).lower()
 if key in cache:return cache[key]
 hh=[];h=hashlib.sha256();n=0
 with path.open() as f:
  for line in f:
   if line.startswith('>'):
    if n:hh.append(h.hexdigest())
    h=hashlib.sha256();n+=1
   else:h.update(line.strip().upper().encode())
 if n:hh.append(h.hexdigest())
 cache[key]={'file_sha256':sha(path),'sequence_sha256':hashlib.sha256('\n'.join(sorted(hh)).encode()).hexdigest()}
 return cache[key]
for i,c in enumerate(annotation):
 if c['scaffolds_file']:continue
 protein=Path(c['protein_file']);genome=protein.with_name(protein.name.replace('.proteins.fa','.fasta').replace('_proteins.fa','.fasta'))
 if genome.exists():
  h=seqhash(genome);c['scaffolds_file']=str(genome);c['scaffolds_file_sha256']=h['file_sha256']
  c['exact_canonical_sequence_matches']=[s for s in c['samples_from_path'] if h['sequence_sha256']==canonical[s]]
 if (i+1)%25==0:print('Annotation assembly links checked',i+1,flush=True)
save(OUT/'data/annotation_candidates.json',annotation)
save(OUT/'data/current_input_hashes.json',cache)
for m in manifest:m['exact_assembly_annotation_candidates']=sum(m['sample'] in a['exact_canonical_sequence_matches'] for a in annotation)
save(OUT/'data/analysis_eligibility.json',manifest)
summary['labels_with_exact_assembly_annotation_candidates']=sum(m['exact_assembly_annotation_candidates']>0 for m in manifest)
save(OUT/'reports/preparation_validation.json',summary)
relfile=OUT/'data/annotation_assembly_relationships.json'
relations=load(relfile) if relfile.exists() else []
byrel={(r['sample'],r['protein_file']):r['relation'] for r in relations}
for m in manifest:
 cc=[a for a in annotation if (m['sample'] in a['exact_canonical_sequence_matches'] or byrel.get((m['sample'],a['protein_file']))=='unchanged_contig_subset') and not a['counts'].get('duplicate_ids',0) and not a['counts'].get('non_amino_acid_records',0) and a['unmapped_proteins']==0]
 cc.sort(key=lambda a:(m['sample'] in a['exact_canonical_sequence_matches'],a['stage']=='annotate',a['mtime_ns']),reverse=True)
 proposals.append({'sample':m['sample'],'preferred_candidate':cc[0]['protein_file'] if cc else None,'assembly_relationship':byrel.get((m['sample'],cc[0]['protein_file'])) if cc else None,'alternative_candidates':[c['protein_file'] for c in cc[1:]],'selection_rule':'exact or unchanged-subset contigs; valid IDs/amino-acid alphabet/GFF parents; prefer exact, annotated, then newest timestamp','status':'proposal only; excluded-contig effects, CDS/translation validation and annotation provenance remain open' if cc else 'manual assembly/GFF link resolution or new annotation needed'})
save(OUT/'data/annotation_release_proposals.json',proposals)
save(OUT/'data/annotation_release_exceptions.json',[p for p in proposals if p['preferred_candidate'] is None])
tools={'checked':'2026-11-18','busco':{'version':'6.1.0','python':'/home/benjamin/.local/share/Achromobacter_ANI/envs/busco/bin/python3.12','executable':'/home/benjamin/.local/share/Achromobacter_ANI/envs/busco/bin/busco','version_and_help_verified':True},'iqtree':{'version':'2.4.0','executable':'/home/benjamin/miniforge3/envs/klebsiella_recombination/bin/iqtree2','version_and_help_verified':True,'original_core_version':'3.1.1','exact_version_reproduction':'requires locating/pinning 3.1.1'},'signalp':{'local_version_and_models':'not verified; required before execution','CLI_documentation':'https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md'},'resource_launcher':{'path':'/usr/bin/systemd-run','cgroup_v2_present':True,'user_manager_report':'degraded','shell_syntax_pass':True,'effective_scope_limits':'must be verified inside each job; not claimed runtime-tested'}}
save(OUT/'data/tool_preflight.json',tools)
actions=load(OUT/'data/outstanding_gate_actions.json')
lines=['TERREUS G0-G4 PREPARATION STATUS — 2026-11-18','','Preparation complete for the available evidence. Scientific gates remain open','where reruns or provenance/metadata decisions are required. No new inference ran.','',f'Canonical assemblies: {len(manifest)} labels / 72 sequence groups; current hashes verified.',f'Core alignment: 4,661 exact partitions, 9,012,257 columns, 66 tips.',f'BUSCO summaries inventoried: {summary["busco_summaries"]}; parse errors: {summary["busco_parse_errors"]}.',f'Labels with provisional matching Eurotiales genome results: {summary["labels_with_provisional_harmonised_genome_busco"]}/77.',f'Unique missing genome BUSCO jobs prepared: {len(queue)}.',f'Labels with candidate raw-read files: {summary["labels_with_candidate_raw_reads"]}/77 (pair/provenance verification pending).',f'Protein annotation candidates audited: {len(annotation)}.',f'Labels with proposed releases after basic checks: {sum(p["preferred_candidate"] is not None for p in proposals)}/77.',f'Incomplete gate actions covered: {len(actions)}.',f'Package validation: {len(checks)} checks passed.','', 'Prepared files:', 'RUNBOOK.txt — exact closure criteria, sequence of work and remaining decisions.', 'data/analysis_eligibility.json — per-label assembly, reads, annotation and BUSCO status.', 'data/annotation_release_proposals.json — candidates, not silently approved versions.', 'data/outstanding_gate_actions.json — each incomplete prior gate and closure evidence.', 'data/busco_missing_genome_queue.json — argument arrays, hashes and prerequisites.', 'data/busco_run_registry.json — raw parameters, lineage, engine and input evidence.', 'data/busco_lineage_lock.json — all local Eurotiales dataset file hashes.', 'inputs/core4661_partitions.nex — recovered non-lexical locus order.', 'scripts/export_primary_proteins.py — GFF-parent selection; retains paralogous genes.', 'scripts/run_one_busco.py — dry run by default, verifies hashes/limits on execution.', 'scripts/run_one_bounded.sh — future resource-window launcher; not run here.', '', 'Unresolved dependencies:', '- Biological aliases, phenotype field definitions and missing assembly/downsampling records.', '- Historical BUSCO input identity; current-file matching alone remains provisional.', '- Assembly-coordinate/CDS translation validation before freezing annotation releases.', '- Exact original IQ-TREE 3.1.1 and downstream predictor/model/database locks.', '- Heavy support/sensitivity, corrected secretome/SNP and functional runs.', '', 'No Panaroo job or global memory setting was changed. Prior outputs remain intact.']
(OUT/'STATUS.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
files=[p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
save(OUT/'manifest.json',{'files':[{'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files],'no_analyses_launched':True})
archive=OUT.parent/'terreus_G0-G4_preparation_2026-11-18_v1.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for p in OUT.rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(OUT))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print('\n'.join(lines[:15]));print('ZIP',archive,archive.stat().st_size)
