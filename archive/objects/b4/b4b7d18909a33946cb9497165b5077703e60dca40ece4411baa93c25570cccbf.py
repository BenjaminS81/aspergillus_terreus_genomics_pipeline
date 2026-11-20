from pathlib import Path
import json, csv, hashlib, shutil
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1]
sources=[]
def read(rel):
    p=ROOT/rel
    ident=f'S{len(sources)+1:02d}'
    dst=OUT/'data'/'sources'/ident/p.name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
    sources.append({'id':ident,'source':str(p),'snapshot':str(dst.relative_to(OUT)).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    return json.loads(p.read_text(encoding='utf-8'))
audit=read('review_reports/2026-11-10/global_status_and_gates.json')
status=read('comparative_visualization/2026-11-11/status_summary.json')
validation=read('comparative_visualization/2026-11-11/validation_summary.json')
rank=read('comparative_visualization/2026-11-11/core66_combinations/ranking.json')
boot=read('comparative_visualization/2026-11-11/core66_combinations/bootstrap_RF_summary.json')
ordn=read('comparative_visualization/2026-11-11/core66_combinations/ordination.json')
metrics=read('comparative_visualization/2026-11-13/data/metric_rows.json')
palette=read('comparative_visualization/2026-11-11/analysis/taxon_palette.json')
ann=read('comparative_visualization/2026-11-11/analysis/core66_tree_annotations.json')
marker_ann=read('comparative_visualization/2026-11-13/marker_provenance/data/provenance_annotations.json')
clades=read('comparative_visualization/2026-11-13/marker_provenance/data/type_clade_key.json')
summary=read('comparative_visualization/2026-11-13/marker_provenance/data/tree_summary.json')
plainval=read('comparative_visualization/2026-11-13/plain_views/validation.json')
recovery=read('taxonomy_reference_build/2026-11-10/status_summary.json')
read('review_reports/2026-11-10/concatenation_order_verification.json')
read('comparative_visualization/2026-11-11/analysis/genome_type_anchor_audit.json')
lov=read('comparative_visualization/2026-11-11/lovastatin/visualization_data.json')
core={a['tip']:a for a in ann if a['tree']=='core_genome'}
config=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\config\samples.tsv')
shutil.copy2(config,OUT/'data'/'project_samples.tsv')
project={x['sample_id']:x for x in csv.DictReader(config.open(),delimiter='\t')}
sources.append({'id':'S18','source':str(config),'snapshot':'data/project_samples.tsv','sha256':hashlib.sha256(config.read_bytes()).hexdigest()})
gates=[]
def gate(family,run,quality,done,remaining,use,src,gate_id):
    gates.append({'analysis':family,'execution':run,'quality':quality,'completed_evidence':done,'remaining_gate':remaining,'stakeholder_use':use,'evidence':src,'gate':gate_id})
gate('Cohort / versions','Complete manifest','Conditional','77 assembly labels; 72 distinct canonical assemblies; 66 core-tree tips. Hash-based sequence metrics are available.','Reconcile aliases to biological cultures; freeze annotation dependencies and analysis-specific eligibility.','Use label counts; do not claim 77 independent isolates.','S02, S03, S07','G0')
gate('Read QC / assembly','Partial','Open','35 selected fastp JSONs and 33 scaffold FASTAs in original audit.','Match reads to exact assemblies; preserve assembly logs and missing downsampling seed/fraction. Read-only QC is N/A if no reads.','Retain audited metrics only; omit unverified read panels.','S01','G1')
gate('Genome / transcript BUSCO','Partial','Open','66 genome and 59 transcript rows in expanded QC; 65 separate Aspergillus BUSCO summaries. Only 60 exact-name core overlaps.','Reconcile assembly hashes, aliases, mode and lineage; run genuinely missing eligible QC. Transcript and genome BUSCO are separate.','Do not present a single harmonised all-genome BUSCO panel yet.','S01','G1')
gate('Full-cohort sequence metrics','Complete','Pass','77 labels / 72 sequences: length, contigs, N50, GC and N metrics; all-genome Circos refreshed Nov 13.','Preserve assembly hashes; this does not close BUSCO or read-quality gates.','All-genome Circos and workbook coverage table.','S07, S13','G5')
gate('Annotation release','Partial','Open','72 annotate summaries and 74 predict summaries were found; these are file counts.','Select one version per assembly; validate primary protein per gene, fragmentation and model integrity.','Annotation-dependent claims remain conditional.','S01','G0 / G4')
gate('Orthogroups','Partial','Conditional','14,964 classic OGs; 6,090 universal; 4,661 strict single-copy; 8,003 in at least 95%. Membership audit passed.','Recover exact protein inputs and completion provenance; finish species-tree/HOG stages if needed.','Descriptive occupancy only; not a completed HOG/species-tree analysis.','S01','G4')
gate('Core concatenation','Complete','Pass','4,661 blocks x 66 samples exactly reconstruct all 9,012,257 columns in recovered locus order.','Preserve recovered coordinates; validate partitions before partitioned reanalysis.','Alignment integrity passed; orthology assumptions remain separate.','S01, S15','G3')
gate('Core-genome ML tree','Complete','Exploratory','IQ-TREE 3.1.1; GTR+F+R10; 66 tips; original October 17 topology retained.','Add bootstrap / SH-aLRT, rooting checks and model/alignment sensitivity; investigate 10 composition failures.','Klebsiella-style figure explicitly labelled exploratory.','S01, S02','G3')
gate('Correct BenA / CaM / RPB2 recovery','Complete','Conditional','198 / 198 core-genome locus calls recovered. Legacy BenA paralog was replaced; IFO6365 BenA is partial.','Keep locus coordinates, all-hit adjudication and partial-locus flags; absence of strong duplicate hits does not prove single copy.','Use corrected November alignments only.','S14, S03','G3')
gate('Seven marker combinations','Complete','Conditional','All singles, pairs and triple: supported reference-inclusive, de-novo core66 and reference-pruned views.','Retain partition/window provenance and unresolved taxon boundaries; verify conflicts before formal species assignments.','Current three-gene tree and 66-tip comparisons.','S03, S12','G3')
gate('RF / NMDS comparison','Complete','Conditional','7,000 bootstrap RF comparisons; identical 66 tips. Triple nRF=0.635; NMDS stress-1=0.0672.','Support the core backbone; test rooting/model/sample sensitivity before treating proximity as taxonomic accuracy.','Comparison figure; intervals are bootstrap percentiles versus a fixed core tree.','S04, S05, S06','G3')
gate('Type / source annotation','Complete','Conditional','283 combined-tree tips: 76 genome-derived and 207 reference; 31 reference type entries plus 3 genome type anchors.','Keep six conflicting genome type claims withheld; distinguish repeated entries from unique cultures and mixed-type clades.','Gold stars and separate provenance/type rings; no forced sensu-stricto names.','S10, S11, S12, S16','G3 / G5')
gate('Lovastatin locus architecture','Complete screen','Conditional','77 labels / 72 assemblies; 101 pieces. 45 nine-gene, 11 partial/divergent, 21 without convincing co-localised block.','Validate CDS integrity, fragmented loci, copy number and read depth where available; experimental function/production remains unknown.','Stacked homolog architecture with flanks; missing hits are not confirmed absence.','S02, S03, S17','G4')
gate('SignalP / secretome','Legacy run exists','Hold','67 files; 26,950 inspected cleavage motifs contain only A/C/G/T; saved recipe uses transcripts.','Rerun validated amino-acid inputs with hashes, model/version and output checks; replace dependent summaries.','Exclude legacy secretome figures from stakeholder set.','S01','G2')
gate('EffectorP','Legacy run exists','Conditional','67 whole-proteome outputs; retained identifiers passed inspected mapping check.','Integrate corrected secretome and transmembrane filtering; distinguish candidates from validated function.','No claims of secreted effectors or human virulence.','S01','G4')
gate('dbCAN','Legacy run exists','Conditional','67 protein-mode overview files.','Separate genes, domains and families; parse multi-domain calls; pin database and evidence threshold.','Defer comparative family heatmap until counting policy is fixed.','S01','G4')
gate('antiSMASH / BiG-SCAPE','Complete clustering','Conditional','BiG-SCAPE 2.0.3: 5,321 records, glocal cutoff 0.3; matrix 67 x 239 BGC families.','Canonical records/IDs, edge-fragment flags and annotation versions; cutoff sensitivity.','BGC families are not genome-wide orthogroups.','S01','G4')
gate('Functional clustering code','Legacy outputs exist','Hold','Stored PCA/dendrogram figures exist; stale OG mapping and saved KeyError were found.','Clean reproducible execution; correct distances/Newick conversion and justified scaling/clusters.','Exclude old functional dendrograms from release.','S01','G4 / G5')
gate('Legacy SNP tree','Complete legacy run','Hold','44 samples; 4,062,390 variable columns; GTR+ASC; 1,000 UFBoot. Missing calls were recoded to reference.','Rebuild callability-aware genotypes/masks; inspect reference bias and within-species eligibility.','Exclude old SNP distances/tree despite completed bootstrap.','S01','G2')
gate('SV / CNV','Partial','Open','43 SV finished, 44 CNV finals, 43 integrated/annotated outputs in audit; cohorts differ.','Read-eligible manifest, stage exit checks, haploidy, depth/repeat assessment and caller concordance.','Read-supported cases only; no all-genome absence claims.','S01','G1 / G4')
gate('Minimum spanning tree','Graph validated','Conditional','67 nodes / 66 edges; minimum stored-matrix weight independently verified as 4,127,821.','Resolve source distances, callable denominator and cohort before redrawing.','Exclude until distance provenance closes; no transmission inference.','S01','G0 / G2 / G5')
gate('Metadata / phenotypes','Partial','Open','12 environmental isolates: MG1=3, MG2=3, MG3=5, unknown=1; R22/R23/R24 metadata incomplete.','Validate metadata joins, assay units/method/endpoints and -1/-2 columns; retain unknowns.','MG annotations only; no resistance categories or causal associations.','S01, S18','G0 / G5')
gate('Ancestral-state / evo-devo','Not run','Not run','Five approaches discussed; no fitMk, stochastic-mapping or phytools result exists in this release.','Supported rooted backbone, validated characters/absence, alias handling and adequate independent transitions.','Next-analysis plan only; no ancestral-state figure.','Discussion; S02','G3 / G4')
next_rows=[
 ['Ancestral block state','Not run','fitMk ER versus ARD; stochastic character mapping','Confirmed absent / partial / complete / unknown states; supported rooted tree','Posterior state pies and transition uncertainty'],
 ['Architecture evolution','Not run','Ordered multistate Mk; optional hidden-rate model','Reliable gene order, strand, contig-break and model integrity calls','Rearrangement/state map with uncertain assemblies flagged'],
 ['Module co-evolution','Not run','Pagel dependent versus independent discrete models','Enough independent transitions; no redundant gene-pair tests','Module dependency estimates; not regulatory causality'],
 ['Cluster / species discordance','Not run','Gene-tree reconciliation; duplication/loss versus transfer models','Ortholog identity, appropriate roots, gene-tree support and contamination checks','Reconciliation / synteny evidence; RF alone does not prove transfer'],
 ['Lineage exception atlas','Not run','Phylogeny-ordered state and architecture matrix','Culture-level de-duplication; technical missingness separated from biological variation','Prioritised validation candidates; no clinical-risk score']]
genomes=[]
for m in metrics:
    a=core.get('LOCAL_'+m['sample']);p=project.get(m['sample'],{})
    genomes.append({**m,'core_tree':bool(a),'reviewed_species':a['species_ring'] if a else 'Not in core tree','project_group':'Environmental study' if p and m['sample'] not in ['R22','R23','R24'] else 'Project R22-R24' if p else 'Other genome','major_genotype':p.get('major_genotype','unknown').replace('none','unknown'),'type_culture':a.get('genome_type_culture','') if a else ''})
figures=[
 ['F01','Three-gene reference tree','01_three_gene_tree.pdf','283 tips; circular + vertical + complete keys','Sequence provenance, type anchors and taxonomic context','Conditional: unresolved species and mixed-type clades retained'],
 ['F02','Core-genome tree','02_core_genome_tree.pdf','66 tips; circular + uncompressed linear + legend','Klebsiella-style metadata rings and study markers','Exploratory: no bootstrap / SH-aLRT for core topology'],
 ['F03','Lovastatin architecture','03_lovastatin_linear.pdf','77 labels / 101 pieces; up to 20 kb flanks','All homolog windows stacked with coordinate annotations','Homology does not demonstrate intact biosynthesis or production'],
 ['F04','Whole-genome Circos','04_whole_genome_circos.pdf','77 labels / 72 canonical sequences','Assembly-wide composition and contiguity','Display concatenates length-sorted contigs; not chromosome synteny'],
 ['F05','Marker choice comparison','05_marker_comparison.pdf','Seven choices on exact same 66 tips','RF ranking + bootstrap percentiles + NMDS','Similarity to an unsupported backbone is not taxonomic validation']]
package={'release':'2026-11-17_v1','review_date':'2026-11-17','evidence_window':'Original audit Nov 10; analytical corrections Nov 10-11; latest figure exports Nov 13-14. No biological reruns in this release.','gates':gates,'sources':sources,'ranking':rank,'bootstrap':boot,'ordination':ordn,'genomes':genomes,'palette':palette,'core_annotations':list(core.values()),'marker_annotations':[a for a in marker_ann if a['tree']=='combined_all_strains'],'type_clades':[c for c in clades if c['tree']=='combined_all_strains'],'next_analyses':next_rows,'figures':figures,'project_metadata':project}
(OUT/'data'/'package.json').write_text(json.dumps(package,indent=2),encoding='utf-8')
(OUT/'data'/'source_manifest.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
for name,rows in [('analysis_gates',gates),('genome_coverage',genomes),('marker_ranking',rank)]:
 with (OUT/'data'/(name+'.tsv')).open('w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('Prepared',len(gates),'gate rows;',len(genomes),'genomes;',len(sources),'sources')
