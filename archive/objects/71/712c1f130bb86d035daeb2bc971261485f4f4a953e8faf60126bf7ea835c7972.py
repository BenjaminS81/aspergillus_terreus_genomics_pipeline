"""Read-only source audit and dependency preparation. Does not launch analyses."""
from pathlib import Path
import sys, os, json, hashlib, subprocess, re, collections, datetime, csv, shutil
ROOT=Path(__file__).resolve().parents[3]; OUT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
PROJECT=Path('D:/Illumina/David_Eisele/terreus_genotyping')
for folder in ['data','reports','scripts','inputs','logs']:(OUT/folder).mkdir(exist_ok=True)
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,v):(OUT/p).write_text(json.dumps(v,indent=2,ensure_ascii=False),encoding='utf-8')
def wpath(p):
 p=str(p)
 if re.match(r'^/mnt/[a-z]/',p):return Path(p[5].upper()+':/'+p[7:])
 return Path(p)
def lpath(p):
 p=Path(p).absolute().as_posix();return '/mnt/'+p[0].lower()+p[2:]
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def fasta(p):
 name=None;parts=[]
 with Path(p).open(encoding='utf-8-sig') as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,''.join(parts).upper()
    name=line[1:].strip();parts=[]
   else:parts.append(line.strip())
 if name is not None:yield name,''.join(parts).upper()
cache=load(OUT/'data/current_input_hashes.json') if (OUT/'data/current_input_hashes.json').exists() else {}
def seqhash(p):
 p=Path(p);key=str(p.resolve()).lower()
 if key in cache and (not p.is_file() or cache[key].get('bytes')!=p.stat().st_size or cache[key].get('mtime_ns')!=p.stat().st_mtime_ns):cache.pop(key)
 if key not in cache:
  if not p.is_file():return None
  hashes=[];bases=0;ids=set();dups=0
  for h,s in fasta(p):
   k=h.split()[0];dups+=k in ids;ids.add(k);bases+=len(s);hashes.append(hashlib.sha256(s.encode()).hexdigest())
  cache[key]={'file_sha256':sha(p),'sequence_sha256':hashlib.sha256('\n'.join(sorted(hashes)).encode()).hexdigest(),'sequences':len(hashes),'bases':bases,'duplicate_ids':dups,'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}
 return cache[key]

source=load(ROOT/'taxonomy_reference_build/2026-11-10/analysis/genome_sources.json')
groups=load(ROOT/'comparative_visualization/2026-11-13/data/genome_manifest.json')
metrics=load(ROOT/'comparative_visualization/2026-11-13/data/genome_metrics.json')
bylabel={s:g for g in groups for s in g['samples']}
assert len(source)==77 and len(groups)==72 and set(bylabel)=={x['sample'] for x in source}
paths=[]
# Bound discovery to final outputs and raw-read directories; do not walk multi-million
# file intermediate BLAST/annotation trees or unrelated jobs.
scan_roots=[PROJECT/'terreus_qc_check',PROJECT/'results/0_2_2_busco',PROJECT/'results/0_2_repeat_mask_funannotate',PROJECT/'results1/0_2_repeat_mask_funannotate',PROJECT/'results2/0_2_repeat_mask_funannotate',PROJECT/'data/raw']
for scan_root in scan_roots:
 for dirname,dirs,files in os.walk(scan_root):
  dirs[:]=[d for d in dirs if not (d.endswith('_misc') or d.startswith('run_') or d in {'busco_downloads','vendor','tmp','work','interproscan','databases','logs','read_mapping','mapping'})]
  for name in files:
   if name.startswith('short_summary.specific.') or name.endswith(('.proteins.fa','_proteins.fa','.transcripts.fa','_transcripts.fa','.summary.json','.fastq.gz','.fq.gz','.fastq','.fq')):paths.append(Path(dirname)/name)
save('data/relevant_file_inventory.json',[str(p) for p in paths if any(x in p.name for x in ['short_summary.specific.','.proteins.fa','.transcripts.fa','_proteins.fa','_transcripts.fa','.summary.json'])])
print('Inventoried',len(paths),'paths; verifying 72 assembly groups',flush=True)
canonical=[]
for i,g in enumerate(groups):
 p=wpath(g['subject']);facts=seqhash(p)
 assert facts and facts['sequence_sha256']==g['sequence_sha256'],f'Canonical assembly changed: {p}'
 canonical.append({**g,'windows_path':str(p),'current_file_sha256':facts['file_sha256'],'hash_verified':True,'biological_culture_equivalence':'not established by sequence identity'})
 if (i+1)%20==0:print('Assembly hashes verified',i+1,flush=True)
save('data/canonical_assemblies.json',canonical)
save('data/sequence_identical_label_groups.json',[g for g in canonical if len(g['samples'])>1])

# Preserve every BUSCO run separately. Hash current source content, never infer historical provenance.
buscos=[]
for p in paths:
 if not (p.name.startswith('short_summary.specific.') and p.suffix=='.json'):continue
 try:
  d=load(p);pa=d.get('parameters',{});lin=d.get('lineage_dataset',{});ip=wpath(pa.get('in',''))
  mode=pa.get('mode','');mode='transcriptome' if 'tran' in mode else 'genome' if 'genome' in mode else 'proteins' if 'prot' in mode else mode
  facts=seqhash(ip) if ip.is_file() else None
  buscos.append({'summary':str(p),'summary_sha256':sha(p),'input':str(ip),'input_current':facts,'mode':mode,'lineage':lin,'versions':d.get('versions',{}),'parameters':pa,'results':d.get('results',{}),'historical_input_hash':'not retained in summary; current content match alone is provisional'})
 except Exception as e:buscos.append({'summary':str(p),'parse_error':str(e)})
save('data/busco_run_registry.json',buscos)
lineage=PROJECT/'terreus_qc_check/busco_downloads/lineages/eurotiales_odb12.2'
assert (lineage/'dataset.cfg').is_file()
if not (OUT/'data/busco_lineage_lock.json').exists():
 lfiles=sorted(p for p in lineage.rglob('*') if p.is_file())
 save('data/busco_lineage_lock.json',{'path':str(lineage),'files':[{'path':str(p.relative_to(lineage)),'sha256':sha(p)} for p in lfiles]})
print('BUSCO summaries reconciled',len(buscos),flush=True)

raw=[p for p in paths if p.name.endswith(('.fastq.gz','.fq.gz','.fastq','.fq')) and 'raw' in p.parts]
proteins=[p for p in paths if p.name.endswith(('.proteins.fa','_proteins.fa')) and any(x in p.parts for x in ['annotate_results','predict_results'])]
transcripts=[p for p in paths if p.name.endswith(('.transcripts.fa','_transcripts.fa')) and any(x in p.parts for x in ['annotate_results','predict_results'])]
alias_candidates={s:{s} for s in bylabel}
for composite in bylabel:
 if re.search(r'_(?:SRR|ERR|DRR)\d+$',composite):
  culture,accession=composite.rsplit('_',1)
  alias_candidates[composite].update([culture,accession])
  if accession in alias_candidates:alias_candidates[accession].add(culture)
  if culture in alias_candidates:alias_candidates[culture].add(accession)
def belongs(p,s):return any(k in p.parts or re.match(re.escape(k)+r'(?:_|\.|$)',p.name) is not None for k in alias_candidates[s])
save('data/path_alias_candidates.json',{'policy':'path-matching candidates derived from explicit composite labels; not proof of biological identity; require sequence/GFF validation','aliases':{k:sorted(v) for k,v in alias_candidates.items()}})
annotation=[]
for p in proteins:
 samples=[s for s in bylabel if belongs(p,s)]
 if not samples:continue
 counts=collections.Counter();ids=set();genes=set();unmapped=0
 gff=p.with_name(p.name.replace('.proteins.fa','.gff3').replace('_proteins.fa','_genes.gff3'))
 mapping={}
 if gff.exists():
  for line in gff.read_text().splitlines():
   if line.startswith('#'):continue
   fields=line.split('\t')
   if len(fields)!=9 or fields[2] not in ['mRNA','transcript']:continue
   a=dict(x.split('=',1) for x in fields[8].split(';') if '=' in x)
   if 'ID' in a and 'Parent' in a:mapping[a['ID']]=a['Parent']
 for h,s in fasta(p):
  pid=h.split()[0];counts['proteins']+=1;counts['residues']+=len(s);counts['duplicate_ids']+=pid in ids;ids.add(pid)
  counts['internal_stop']+='*' in s.rstrip('*');counts['terminal_stop']+=s.endswith('*')
  counts['non_amino_acid_records']+=bool(set(s)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO*'))
  counts['nucleotide_only_records']+=bool(s) and set(s)<=set('ACGTN');counts['not_starting_M']+=not s.startswith('M')
  if pid in mapping:genes.add(mapping[pid])
  else:unmapped+=1
 gfasta=p.with_name(p.name.replace('.proteins.fa','.scaffolds.fa').replace('_proteins.fa','_scaffolds.fa'))
 if not gfasta.exists():gfasta=p.with_name(p.name.replace('.proteins.fa','.fasta').replace('_proteins.fa','.fasta'))
 genome=seqhash(gfasta) if gfasta.exists() else None
 matches=[s for s in samples if genome and genome['sequence_sha256']==bylabel[s]['sequence_sha256']]
 annotation.append({'protein_file':str(p),'samples_from_path':samples,'file_sha256':sha(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns,'stage':'annotate' if 'annotate_results' in p.parts else 'predict','gff':str(gff) if gff.exists() else None,'gff_sha256':sha(gff) if gff.exists() else None,'scaffolds_file':str(gfasta) if gfasta.exists() else None,'exact_canonical_sequence_matches':matches,'counts':dict(counts),'mapped_genes':len(genes),'unmapped_proteins':unmapped,'primary_isoform_policy':'one longest valid protein per explicit GFF parent; deterministic ID tie-break; preserve discarded mapping; do not collapse paralogous genes','release_status':'candidate; requires CDS/translation and assembly-coordinate checks'})
save('data/annotation_candidates.json',annotation)

manifest=[];queue=[];buscoenv='/home/benjamin/.local/share/Achromobacter_ANI/envs/busco'
for item in source:
 s=item['sample'];g=bylabel[s]
 label_facts=seqhash(wpath(item['subject']))
 assert label_facts and label_facts['sequence_sha256']==g['sequence_sha256'],f'Assembly alias changed: {s}'
 candidates=[a for a in annotation if s in a['samples_from_path']]
 reads=[p for p in raw if belongs(p,s)]
 trs=[p for p in transcripts if belongs(p,s)]
 matches=[b for b in buscos if b.get('mode')=='genome' and b.get('input_current') and b['input_current']['sequence_sha256']==g['sequence_sha256'] and b.get('lineage',{}).get('name')=='eurotiales_odb12.2' and b.get('results',{}).get('n_markers')==4240]
 row={'sample':s,'representative':g['representative'],'sequence_sha256':g['sequence_sha256'],'assembly':str(wpath(item['subject'])),'in_core_tree':item['in_core_tree'],'raw_read_candidates':[str(p) for p in reads],'read_eligibility':'candidate reads; pair/lane/assembly provenance not yet validated' if reads else 'no reads found in project raw directories; read-only analyses N/A pending alias review','protein_candidates':len(candidates),'transcript_candidates':[str(p) for p in trs],'exact_assembly_annotation_candidates':sum(s in a['exact_canonical_sequence_matches'] for a in candidates),'busco_genome_matching_current_input':[b['summary'] for b in matches],'busco_status':'provisional existing same-sequence/same-lineage result; historical hash unresolved' if matches else 'missing harmonised same-sequence genome result','metadata_status':'join and biological-culture validation pending'}
 row['assembly_file_sha256']=label_facts['file_sha256']
 manifest.append(row)
 if not matches and s==g['representative']:
  queue.append({'job':'BUSCO_genome_'+s,'status':'prepared_not_run','requires':['canonical hash verified','resource wrapper','tool and lineage lock match'],'input_windows':str(wpath(g['subject'])),'expected_file_sha256':seqhash(wpath(g['subject']))['file_sha256'],'argv':[buscoenv+'/bin/python3.12',buscoenv+'/bin/busco','--in',g['subject'],'--out','genome_'+s,'--out_path',lpath(OUT/'runs/busco'),'--mode','genome','--lineage_dataset',lpath(lineage),'--offline','--cpu','1','--miniprot','--opt-out-run-stats'],'environment':{'PATH':buscoenv+'/bin:/usr/local/bin:/usr/bin:/bin'}})
save('data/analysis_eligibility.json',manifest);save('data/busco_missing_genome_queue.json',queue)
save('data/transcript_busco_candidates.json',[{'path':str(p),'samples':[s for s in bylabel if belongs(p,s)],'status':'select annotation release first; never substitute for genome mode'} for p in transcripts if any(belongs(p,s) for s in bylabel)])

# Recover the exact non-lexical partition order; do not silently sort loci.
partition=load(ROOT/'review_reports/2026-11-10/concatenation_order_verification.json');parts=partition['parts']
assert partition['mismatch_count']==0 and partition['continuous_complete_coverage']
assert len(parts)==4661 and parts[0]['start_1based']==1 and parts[-1]['end_inclusive']==9012257
assert all(a['end_inclusive']+1==b['start_1based'] for a,b in zip(parts,parts[1:]))
(OUT/'inputs/core4661_partitions.nex').write_text('#nexus\nbegin sets;\n'+''.join(f"  charset {p['orthogroup']} = {p['start_1based']}-{p['end_inclusive']};\n" for p in parts)+'end;\n')
save('inputs/core4661_coordinates.json',partition)

# Frozen evidence and one action record for every incomplete gate through G4.
pkg=load(ROOT/'stakeholder_release/2026-11-17_v1/data/package.json')
actions=[]
for i,g in enumerate(pkg['gates']):
 if g['quality']=='Pass' or not any(x in g['gate'] for x in ['G0','G1','G2','G3','G4']):continue
 actions.append({'id':f'A{i+1:02d}',**g,'preparation_status':'specified; analytical gate remains '+g['quality'],'closure_evidence_required':g['remaining_gate'],'execution_status':'not launched'})
save('data/outstanding_gate_actions.json',actions)
for name in ['analysis_gates.tsv','source_manifest.json']:
 shutil.copy2(ROOT/'stakeholder_release/2026-11-17_v1/data'/name,OUT/'data'/('previous_'+name))
shutil.copy2(PROJECT/'config/samples.tsv',OUT/'data/project_samples_snapshot.tsv')
save('data/current_input_hashes.json',cache)
summary={'date':datetime.datetime.now().astimezone().isoformat(),'limits':LIMITS,'source_read_only':True,'source_file_paths_indexed':len(paths),'canonical_labels':len(manifest),'canonical_sequence_groups':len(canonical),'all_canonical_sequence_hashes_verified':True,'core_tree_labels':sum(m['in_core_tree'] for m in manifest),'busco_summaries':len(buscos),'busco_parse_errors':sum('parse_error' in b for b in buscos),'labels_with_provisional_harmonised_genome_busco':sum(bool(m['busco_genome_matching_current_input']) for m in manifest),'unique_genome_busco_jobs_prepared':len(queue),'labels_with_candidate_raw_reads':sum(bool(m['raw_read_candidates']) for m in manifest),'protein_release_candidates_audited':len(annotation),'labels_with_exact_assembly_annotation_candidates':sum(m['exact_assembly_annotation_candidates']>0 for m in manifest),'outstanding_action_rows':len(actions),'partition_blocks':len(parts),'partition_columns':partition['length'],'analyses_launched':0,'gate_closure_claimed':False}
save('reports/preparation_validation.json',summary)
print(json.dumps(summary,indent=2),flush=True)
