"""Validate the figures and assemble an offline, searchable review package."""
from circos_review_common import *
LIMITS=constrain_process()
import csv, shutil, html, hashlib, zipfile, collections, re, sys
import pymupdf
D=OUT/'data';F=OUT/'figures';Q=OUT/'qa'
order=readjson(D/'figure_order.json');raw=readjson(D/'visualization_data.json');data={x['sample']:x for x in raw}
metrics=readjson(D/'genome_metrics.json');ann=readjson(D/'sample_annotations.json');inventory=readjson(D/'figure_inventory.json')
tests=[]
def check(name,value):
    tests.append({'check':name,'passed':bool(value)})
    assert value,name
check('77 unique labels in a common order',len(order)==77 and len(set(order))==77 and set(order)==set(data))
check('72 canonical assemblies',len(metrics)==72)
check('all assemblies and region sequences match prior hashes',all(readjson(D/'source_validation.json')[x] for x in ('all_assembly_hashes_match','all_region_hashes_match')))
check('all five figures contain every label',len(inventory)==5 and all(x['sample_count']==77 and set(x['samples'])==set(order) for x in inventory.values()))
check('101 separate segments preserved',len(readjson(D/'region_display_coordinates.json'))==sum(len(x['segments']) for x in raw)==101)
check('whole assembly coordinates cover every base',all(sum(c['length'] for c in m['contig_coordinates'])==m['assembly_bp'] and m['contig_coordinates'][-1]['end_0based_exclusive']==m['assembly_bp'] for m in metrics.values()))
check('all genome bins cover whole assemblies without gaps',all(len(m['bins'])==100 and m['bins'][0]['start']==0 and m['bins'][-1]['end']==m['assembly_bp'] and all(a['end']==b['start'] for a,b in zip(m['bins'],m['bins'][1:])) for m in metrics.values()))
check('contig density accounts for all contigs',all(sum(x['contig_starts'] for x in m['bins'])==m['contigs'] for m in metrics.values()))
check('N50 and composition within valid bounds',all(0<m['N50_bp']<=m['assembly_bp'] and 0<=m['gc_percent']<=100 and 0<=m['n_percent']<=100 for m in metrics.values()))
check('all gene hits lie in their displayed windows',all(p['window_start']<=g['start']<=g['end']<=p['window_end']<=p['contig_length'] for x in raw for p in x['segments'] for g in p['genes']))
check('44 explicitly scoped lovB links',inventory['B_all_genomes_coordinate_circos']['lovB_links']==44)
expected=collections.Counter(x['block_status'] for x in raw)
check('previous evidence categories preserved',expected=={'nine_gene_block':45,'partial_or_divergent_block':11,'no_convincing_colocalized_block':21})
for stem,info in inventory.items():
    doc=pymupdf.open(F/(stem+'.pdf'));check(stem+' is a valid one-page vector PDF',len(doc)==1 and len(doc[0].get_text())>100)
    text=doc[0].get_text()
    check(stem+' contains all long sample identifiers',all(s in text for s in order if len(s)>4))
    check(stem+' has all export formats',all((F/(stem+'.'+ext)).stat().st_size>10000 for ext in ('pdf','svg','png')))
    scale=.4 if 'stacked' in stem else .6
    doc[0].get_pixmap(matrix=pymupdf.Matrix(scale,scale)).save(str(Q/(stem+'_pdf.png')))
    if 'stacked' in stem:
        for label,y in [('top',0),('middle',1800),('bottom',max(0,doc[0].rect.height-450))]:
            doc[0].get_pixmap(matrix=pymupdf.Matrix(1,1),clip=pymupdf.Rect(0,y,doc[0].rect.width,min(y+450,doc[0].rect.height))).save(str(Q/('linear_'+label+'_pdf.png')))
    doc.close()
writejson(D/'validation.json',{'checks':tests,'passed':len(tests),'failed':0,'visual_review':'PDF renders inspected separately before delivery'})
sources=[
 {'title':'Circos: quick-start histograms','url':'https://www.circos.ca/documentation/tutorials/quick_start/histograms/','used':'Position-value tracks, explicit inner/outer radii, fixed scales; applied in all three circular designs.'},
 {'title':'Circos: heat-map tracks','url':'https://www.circos.ca/documentation/tutorials/2d_tracks/heat_maps/lesson','used':'Sequential color tracks for numerical values; used for composition and QC.'},
 {'title':'pyCirclize: genomic examples (E. coli)','url':'https://moshi4.github.io/pyCirclize/circos_plot/','used':'Separate forward/reverse features, GC and GC-skew rings. Adapted to linear fungal contig windows in design A.'},
 {'title':'pyCirclize: comparative genomics examples','url':'https://moshi4.github.io/pyCirclize/comparative_genomics/','used':'Reference/query sectors and links; design B uses only pre-existing validated lovB homolog intervals, not the example whole-genome MUMmer workflow.'},
 {'title':'pyCirclize: Track API','url':'https://moshi4.github.io/pyCirclize/api-docs/track/','used':'Checked signatures for rect, arrow, heatmap, text, scatter and plotting. The installed API was also inspected.'}
]
writejson(D/'documentation_and_examples.json',{'review_date':'2026-11-13','sources':sources,'implementation':'pyCirclize 1.10.1; Perl Circos was not invoked','note':'Some direct circos.ca fetches failed; indexed official tutorial content was available. The pyCirclize examples and API pages were directly readable.'})
designs=[
 ('A_lovastatin_regions_circos','A · Lovastatin regions with flanks','Best for local context. Every genome contributes its detected homolog regions, with up to 20 kb flanks, forward/reverse features, query coverage, separate DNA and protein identity, GC and oriented GC skew.','Dense with 77 labels; use the stacked companion for individual gene order. Each contig piece is independent.'),
 ('B_all_genomes_coordinate_circos','B · Whole-assembly Circos','Best for genome-wide context. Sector size follows assembly length. Tracks show GC, skew, Ns, contig-start density, pathway positions and N50. Visible lines connect NIH2624 lovB to 44 other nine-gene-block labels.','Contigs are sorted by length for display. Lines show locus homology, not whole-genome alignment or chromosome adjacency.'),
 ('C_all_genomes_qc_evidence_circos','C · QC and pathway comparison — recommended overview','Best for comparing all samples. Equal sectors show assembly size, GC, N50, contig count, N bases, the nine-gene panel, block status, reported taxon and documented type anchors.','Summary matrix, not genomic coordinates or a phylogenetic tree. Reported names marked * remain unresolved or conflicting.'),
 ('lovastatin_all_genomes_stacked_flanks','Linear · All loci and flanks in one plot','A single vector figure with 101 rows for all 77 labels. Panel genes and reference flank homologs are labelled; native contig intervals, orientation, block status and edge flags accompany each row.','Large-format figure: zoom in or use SVG/PDF. Other annotation-model identifiers are searchable below and in the feature table.'),
 ('genome_key','Genome key and identical-assembly aliases','The same numbering is used in every figure. Identical assembly aliases are retained and documented.','77 labels represent 72 sequence sets; they are not 77 independent biological replicates.')
]
cards=[]
for stem,title,why,caution in designs:
    cards.append(f'<article id="{stem}"><h2>{html.escape(title)}</h2><p>{html.escape(why)}</p><p class="note">{html.escape(caution)}</p><nav>'+''.join(f'<a href="figures/{stem}.{ext}">{ext.upper()}</a>' for ext in ('pdf','svg','png'))+f'</nav><a href="figures/{stem}.svg"><img loading="lazy" src="figures/{stem}.png" alt="{html.escape(title)}"></a></article>')
tablerows=[]
for i,s in enumerate(order):
    d=data[s];m=metrics[d['representative']];a=ann['LOCAL_'+s]
    cells=[f'{i+1:02d}',s,d['taxon'],a['species_ring'],f"{m['assembly_bp']/1e6:.3f}",f"{m['gc_percent']:.3f}",f"{m['N50_bp']:,}",f"{m['contigs']:,}",f"{m['n_percent']:.3f}",d['block_status'].replace('_',' '),str(d['genome_wide_strong_genes'])+'/9',str(d['segment_count']),', '.join(m['samples']) if len(m['samples'])>1 else '',a.get('genome_type_culture','')]
    tablerows.append('<tr>'+''.join('<td>'+html.escape(c)+'</td>' for c in cells)+'</tr>')
feature_rows=[]
with (D/'region_features.tsv').open(encoding='utf-8') as h:
    for r in csv.DictReader(h,delimiter='\t'):
        cells=[r[k] for k in ('sample','segment','contig','feature_kind','gene_or_model','start_1based','end_1based_inclusive','strand','method','identity_percent','query_coverage')]
        feature_rows.append('<tr>'+''.join('<td>'+html.escape(c)+'</td>' for c in cells)+'</tr>')
sourcehtml=''.join(f'<li><a href="{x["url"]}">{html.escape(x["title"])}</a> — {html.escape(x["used"])}</li>' for x in sources)
doc='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Terreus | Circos figure review</title><style>
:root{font-family:system-ui,sans-serif;color:#213c50;background:#f3f6f8}body{margin:0}header,main{max-width:1300px;margin:auto;padding:28px}header{padding-top:45px}h1{font-size:36px;margin-bottom:10px}h2{font-size:24px}p{line-height:1.65;max-width:1050px}.eyebrow{letter-spacing:2px;color:#368e87;font-weight:700}.banner{padding:18px 24px;background:#e3f0ed;border-radius:12px}.cards{display:grid;grid-template-columns:1fr 1fr;gap:24px}article,section{background:white;padding:24px;border-radius:12px;margin-bottom:24px;box-shadow:0 2px 12px #1932420b}article img{width:100%;max-height:850px;object-fit:contain;object-position:top;border:1px solid #e6ecef;margin-top:15px}article:nth-child(4){grid-column:1/-1}article:nth-child(4) img{max-height:1100px}nav{display:flex;gap:12px;flex-wrap:wrap}a{color:#187b7a}nav a{padding:8px 14px;background:#e8f2f2;border-radius:6px;text-decoration:none}.note{color:#617786;font-size:14px}input{padding:12px;width:90%;max-width:600px;font-size:15px;border:1px solid #bbcbd3;border-radius:6px}.scroll{overflow:auto;max-height:650px;margin-top:15px}table{border-collapse:collapse;font-size:13px;width:100%}th{position:sticky;top:0;background:#e8eff3;text-align:left}td,th{padding:9px;border-bottom:1px solid #e2e8eb;white-space:nowrap}tr:nth-child(even){background:#f5f8fa}li{line-height:1.65;margin-bottom:8px}code{background:#f1f4f5;padding:2px 5px}summary{cursor:pointer;font-weight:600}@media(max-width:850px){.cards{display:block}header,main{padding:16px}}
</style><header><div class="eyebrow">ASPERGILLUS · FIGURE REVIEW · 13 NOVEMBER 2026</div><h1>Three circular views, one complete cohort</h1><p>77 assembly labels · 72 distinct assemblies · 101 lovastatin-related contig segments. All sequence hashes match the previous analysis. Existing source data and the November 11 figures were preserved.</p><div class="banner"><b>Suggested pairing:</b> C for an overview of genome quality and pathway evidence, plus the stacked linear figure for precise locus comparison. A adds detailed regional metrics; B shows whole-assembly context.</div></header><main><div class="cards">'''+''.join(cards)+'''</div>
<section><h2>Coverage, evidence and limits</h2><p>The existing screen classifies 45 labels as a contiguous nine-gene block, 11 as partial/split, and 21 as no convincing co-localized block. These correspond to 42, 9 and 21 distinct assemblies. A scattered homolog does not establish a cluster. A missed hit does not establish biological absence, and no plot demonstrates lovastatin production.</p><p>All 77 canonical labels in the current genome/locus analysis are present, including reference and FASTA-only assemblies. This is not every archived assembly variant. Transcript-only files cannot provide genomic flanks or whole-assembly coordinates and are not treated as additional genomes. The full earlier inventory of alternative assemblies was not rerun.</p><p>Names come from the previous taxonomy review: * flags unresolved/conflicting reported labels; [T] marks three documented genome type anchors. No new taxonomic assignments were made. Identical-sequence aliases remain visible for traceability.</p><p>Raw-read depth and BUSCO are omitted from these three drafts because a harmonized, exact-assembly-matched dataset has not been verified across the full cohort. Full-genome gene density, repeat density and orthogroup copy-number tracks are likewise not inferred from inconsistent annotation versions.</p></section>
<section><h2>All-genome metrics</h2><p>Search by sample, taxon, alias or status. <a href="data/genome_metrics.tsv">Download tab-delimited metrics</a>.</p><input id="sampleQuery" placeholder="Search all 77 rows" oninput="filterRows('sampleQuery','sampleTable')"><div class="scroll"><table id="sampleTable"><thead><tr>'''+''.join('<th>'+x+'</th>' for x in ['#','Sample','Reported taxon','Reviewed species','Mb','GC %','N50 bp','Contigs','N %','Block status','Strong genes','Pieces','Identical sequence aliases','Type culture'])+'''</tr></thead><tbody>'''+''.join(tablerows)+'''</tbody></table></div></section>
<section><h2>Exact regional annotation</h2><p>All labelled homologs, reference flank hits and other displayed annotation models. Native coordinates are 1-based inclusive. DNA and protein identities are separate quantities; method is retained for every hit. Genomic intervals are homolog spans, not validated intact coding sequences. <a href="data/region_features.tsv">Download the full feature table with provenance</a>.</p><input id="featureQuery" placeholder="Search sample, gene, contig or model identifier" oninput="filterRows('featureQuery','featureTable')"><div class="scroll"><table id="featureTable"><thead><tr>'''+''.join('<th>'+x+'</th>' for x in ['Sample','Piece','Contig','Evidence','Gene / model','Start','End','Strand','Method','Identity %','Query coverage'])+'''</tr></thead><tbody>'''+''.join(feature_rows)+'''</tbody></table></div></section>
<section><h2>Documentation and examples used</h2><p>The figures were implemented with pyCirclize 1.10.1 and Matplotlib, following Circos track conventions. Perl Circos itself was not invoked. The official examples informed the design; the figures above contain this project's data.</p><ul>'''+sourcehtml+'''</ul></section>
<section><h2>Methods and reproducibility</h2><ul><li>Sequence metrics were recalculated from the 72 canonical FASTAs, using one worker limited to 2 GiB RAM and one CPU core. This worker did not modify the Panaroo job or its settings.</li><li>GC = 100 × (G+C)/(A+C+G+T). N % uses all sequence positions as denominator; the accompanying JSON also retains all-ambiguity %. GC skew = (G-C)/(G+C), undefined only where G+C=0. A reverses the skew sign for reversed display segments; B preserves submitted contig orientation. Neither plot treats skew as a fungal replication-origin inference.</li><li>Regional metrics use non-overlapping 1 kb windows, including a shorter final window. Whole-assembly tracks use 100 near-equal base-count bins per assembly. No gaps are inserted between contigs in B's display coordinates; the contig-start-density track reveals fragmentation.</li><li>N50 is calculated across all FASTA sequence records without a minimum-length cutoff or splitting at Ns. Each assembly's contigs are sorted by decreasing length and then identifier for B; exact offsets are in genome_metrics.json. Artificial 2 kb gaps in A separate extracted fragments and are excluded from sequence metrics.</li><li>Previous BLASTn hits use at least 70% query coverage and 75% nucleotide identity; protein-rescue hits use at least 70% coverage and 45% amino-acid identity. The nine-gene panel is a screening panel, not a complete functional/BGC boundary validation. DNA and amino-acid identity never share a single scale.</li><li>Arrows use genomic spans and orientation. A feature clipped at its directional end has no arrowhead. Matching-coordinate gene models were inherited from the validated November 11 analysis; unverified annotation contexts are not fabricated.</li><li>Full precision data, exact display-coordinate maps, validation results and scripts accompany the figures. Numerical color ranges are fixed; out-of-range values saturate. A uses GC 30-70% and skew +/-0.5; B uses GC 40-60%, skew +/-0.08, N 0-5%, log10(1+contig starts) 0-4 and N50 10^3-10^7 bp.</li></ul><nav><a href="data/validation.json">Validation checks</a><a href="data/source_validation.json">Sequence checks</a><a href="data/figure_inventory.json">Figure coverage</a><a href="REPRODUCE.txt">Reproduction instructions</a></nav></section>
<section><h2>Useful additions after these drafts</h2><table><tr><th>Track / extension</th><th>Required gate</th></tr><tr><td>BUSCO completeness</td><td>Match exact assembly versions, lineage release and analysis mode; leave unavailable values explicitly missing.</td></tr><tr><td>Read depth / locus copy evidence</td><td>Map available raw reads to the matching assembly; normalize against genome depth. Omit samples without reads.</td></tr><tr><td>Full-genome synteny ribbons</td><td>Generate and filter actual genome alignments; select a reference or species subsets to avoid a 77-genome link tangle.</td></tr><tr><td>Repeat / gene / orthogroup density</td><td>Validate annotation coordinates and comparable calling methods for each genome.</td></tr><tr><td>Functional cluster integrity</td><td>Inspect gene models, frameshifts, boundaries and duplicated candidates; production requires separate experimental evidence.</td></tr></table></section></main><script>function filterRows(input,table){const query=document.getElementById(input).value.toLowerCase();document.querySelectorAll('#'+table+' tbody tr').forEach(row=>row.hidden=!row.textContent.toLowerCase().includes(query));}</script></html>'''
(OUT/'circos_review.html').write_text(doc,encoding='utf-8')
instructions='''Terreus Circos review, 2026-11-13

Open circos_review.html. All four primary figures contain all 77 canonical labels.
PDF and SVG are scalable; the linear figure is a single large-format page.

Plot reproduction (from this output folder):
  python -m pip install --target vendor -r requirements.txt
  python scripts/plot_circos_review.py
Optional plot selection: A B C linear key
The renderer uses only the saved data/ inputs. No source FASTAs or raw reads
are required to replot. Exact package versions are in plotting_environment.json.
On Windows the worker enforces its own 2 GiB job limit and one-CPU affinity.

Recalculate metrics only when needed, using prepare_circos_review.py in the
original workspace. This requires the original D: FASTAs and November 11 data.
The recomputation asserts canonical assembly and extracted-window hashes.

The vendor installation and plotting cache are intentionally excluded from
the distributable ZIP. No old result directories were overwritten.
''' 
(OUT/'REPRODUCE.txt').write_text(instructions,encoding='utf-8')
(OUT/'requirements.txt').write_text('pycirclize==1.10.1\nmatplotlib==3.11.2\nnumpy==2.4.6\nbiopython==1.88\npandas==3.1.1\n',encoding='utf-8')
for name in ('circos_review_common.py','prepare_circos_review.py','plot_circos_review.py','finalize_circos_review.py'):
    shutil.copy2(BASE/name,OUT/'scripts'/name)
# Check local links and the table inventories without requiring a web server.
links=re.findall(r'(?:href|src)="([^"]+)"',doc)
check('all local report links resolve',all((OUT/x).exists() for x in links if not x.startswith(('http','#'))))
check('metric table has 77 rows',len(tablerows)==77)
check('feature table is populated',len(feature_rows)>1000)
writejson(D/'validation.json',{'checks':tests,'passed':len(tests),'failed':0,'visual_review':'PDF renders inspected separately before delivery'})
print('VALIDATED',len(tests),'checks;',len(feature_rows),'feature rows; report created',flush=True)
