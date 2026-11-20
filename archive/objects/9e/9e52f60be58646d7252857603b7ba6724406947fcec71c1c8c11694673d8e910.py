"""Evidence-based gate ledger. Missing or merely queued work never becomes PASS."""
from pathlib import Path
import json,collections,datetime
from worker import OUT,PREP,load,save,sha
def read_optional(name,default):
 p=OUT/'reports'/name;return load(p) if p.exists() else default
actions=load(PREP/'data/outstanding_gate_actions.json')
busco=list((OUT/'reports').glob('BUSCO_genome_*.validated.json'))
annotation=read_optional('annotation_release_summary.json',[])
counts=collections.Counter(r['technical_status'] for r in annotation)
functional=read_optional('functional_reanalysis_acceptance.json',{})
core=read_optional('core_support.validated.json',{})
readqc=read_optional('read_qc_stage.json',{})
proteinqc=read_optional('protein_busco_stage.json',{})
comparison=read_optional('new_core_comparison.json',{})
unresolved_annotations=[r['sample'] for r in annotation if r['technical_status']=='ASSEMBLY_ANNOTATION_LINK_UNRESOLVED']
annotation_links=('unresolved links: '+', '.join(unresolved_annotations)) if unresolved_annotations else ('no unresolved links in completed records' if annotation else 'assembly-link validation pending')
updates={
 'A01':('PARTIAL','77 labels / 72 sequence groups frozen; analysis-specific eligibility retained. Biological culture equivalences remain separate from sequence equality.'),
 'A02':('PARTIAL' if readqc else 'OPEN',f"Read QC stage: {readqc if readqc else 'pending'}. Candidate pairing, raw quality and mapping are measured where available; historical downsampling and biological read-to-isolate identity remain separate."),
 'A03':('GENOME_AND_PRIMARY_GENE_SET_PANELS_PASS' if len(busco)==72 and proteinqc else 'PARTIAL',f'{len(busco)}/72 unique genome BUSCO jobs validated (35 missing cases plus 37 historical-provenance verification runs). Primary-protein mode: {proteinqc if proteinqc else "pending"}. Historical transcript-mode runs remain separate from these harmonised panels.'),
 'A05':('PARTIAL',f'CDS/protein validation status: {dict(counts)} across {len(annotation)}/77 labels; {annotation_links}. Subset-contig consequences and release comparability remain explicit.'),
 'A06':('PARTIAL' if (OUT/'reports/orthofinder_input_recovery.json').exists() else 'OPEN','Original OrthoFinder 3.0.1b1 explicitly used DNA mode (-d), not protein mode. Exact DNA inputs are recovered when the recovery report passes; final species-tree/HOG completion and a primary-protein benchmark remain separate.'),
 'A08':('PARTIAL' if core else 'OPEN','Same-version IQ-TREE 3.1.1 support '+('validated; rooting and model/alignment sensitivity still required.' if core else 'queued or awaiting successful capped execution.')),
 'A09':('CONDITIONAL','Corrected 198 core marker calls retained; full locus-coordinate, duplicate/copy and partial-locus adjudication still required for closure.'),
 'A10':('CONDITIONAL','Existing supported marker combinations retained; unresolved taxon boundaries remain.'),
 'A11':('CONDITIONAL',f'Updated-core comparison: {comparison if comparison else "pending; previous RF/NMDS retained"}. Rooting/model/alignment sensitivity acceptance remains separate.'),
 'A12':('CONDITIONAL','Existing supported type/source annotations retained, conflicting genome type claims withheld.'),
 'A13':('CONDITIONAL','Homolog architecture is available; intact CDS, copy number, contig-boundary/absence and read-depth evidence remain to be validated.'),
 'A14':('BLOCKED_DEPENDENCY','Legacy nucleotide-input SignalP results excluded. A usable SignalP executable/model set with appropriate permitted use and a pinned transmembrane predictor have not been located.'),
 'A15':('CONDITIONAL','Legacy whole-proteome counts reconciled where IDs match validated primary proteins; secreted-effector status awaits corrected SignalP/TM evidence.'),
 'A16':('CONDITIONAL','Multi-domain calls split into family/gene/domain units with explicit support. Legacy input/database hashes are unresolved; no claim of a fresh dbCAN inference.'),
 'A17':('CONDITIONAL','BiG-SCAPE record identity, annotation links, edge fragments and cutoff sensitivity remain pending.'),
 'A18':('PARTIAL' if functional else 'OPEN',('Deterministic Jaccard/average-linkage regenerated without a forced cluster count' if functional else 'Functional reanalysis pending')+'; biological use remains conditional on predictor/database provenance.'),
 'A19':('HOLD','Missing-to-reference legacy SNP analysis excluded. Callable masks, validated haploid read calls and within-species reference-bias checks still needed.'),
 'A20':('OPEN','SV/CNV require qualified read/assembly pairs and caller-stage/depth/repeat validation; no all-genome absence claims.'),
 'A21':('OPEN','Legacy MST graph validation is retained; its source distances/callable denominator remain unresolved.'),
 'A22':('DEFERRED_BY_USER','Phenotype units/method/endpoints, -1/-2 definitions, and R22/R23/R24 metadata intentionally skipped. Keep unknown; list required records at handover.'),
 'A23':('DOWNSTREAM_NOT_RUN','ASR awaits a supported rooted backbone, validated character states and independent-culture handling.')}
for a in actions:
 a['current_status'],a['current_evidence']=updates.get(a['id'],('OPEN',a['remaining_gate']))
 a['previous_quality']=a['quality']
save(OUT/'reports/gate_ledger.json',actions)
needed=[
 'Phenotype dictionary: assay name, units, medium, incubation conditions, endpoint/readout, replicate handling, and meanings of -1 and -2.',
 'R22/R23/R24: source/host or environment, collection date and location, biological isolate/culture ID and sequencing-run/assembly link.',
 'Culture/accession equivalences for repeated labels: collection or BioSample records establishing whether sequences represent the same culture.',
 'SRR33623979 assembly/downsampling provenance: original read pair/lane IDs, exact fraction, seed and software command, or an explicit irrecoverable-record designation.',
 'SignalP executable/model bundle and permitted use; transmembrane predictor/model version. Legacy SignalP nucleotide inputs cannot be reused.',
 'Original predictor/orthogroup database versions and input hashes, or replacement runs from the frozen primary proteins.',
]
if unresolved_annotations:needed.append('Assembly-specific annotation models or documented coordinate/sequence reconciliation for: '+', '.join(unresolved_annotations)+'.')
elif not annotation:needed.append('Annotation-link requirements will be refined after the boundary and CDS validation stage completes.')
(OUT/'REMAINING_REQUIREMENTS.txt').write_text('RECORDS AND DEPENDENCIES STILL NEEDED\n\n'+'\n\n'.join(f'{i+1}. {v}' for i,v in enumerate(needed))+'\n')
lines=['TERREUS G0-G4 EXECUTION STATUS','',f'Updated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}',f'New unique-genome BUSCO validations: {len(busco)}/72 (35 missing + 37 provenance verification).',f'Annotation validation: {dict(counts)}.',f'Core support validated: {bool(core)}.','Sequential analysis: BUSCO capped at 20 GB (20,000,000,000 bytes); other stages 4 GiB. One CPU, no swap, 12 GiB headroom guard.','The protected 100,000,000,000-byte service limit is checked, never changed.','Metadata work is deferred by the user. No missing values are imputed.','','Gate actions:']
lines += [f"{a['id']} | {a['gate']} | {a['analysis']} | {a['current_status']}\n  {a['current_evidence']}" for a in actions]
(OUT/'STATUS.txt').write_text('\n'.join(lines)+'\n')
save(OUT/'reports/acceptance_summary.json',{'updated':datetime.datetime.now(datetime.timezone.utc).isoformat(),'actions':len(actions),'status_counts':dict(collections.Counter(a['current_status'] for a in actions)),'all_G0_G4_closed':False,'metadata_deferred_by_user':True})
print('Gate ledger updated; uncompleted, conditional and deferred work remains explicit.')
