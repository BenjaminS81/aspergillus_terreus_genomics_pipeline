from pathlib import Path
import sys,json,hashlib,shutil,zipfile,xml.etree.ElementTree as E,re,html,collections
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
import pymupdf
P=json.loads((OUT/'data/package.json').read_text());F=OUT/'figures';Q=OUT/'qa'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name in ['circos_review_common.py','terreus_tree_helpers.py','plot_plain_terreus_views.py']:
 shutil.copy2(ROOT/name,OUT/'scripts'/name)
for name in ['core_genome_original.treefile','combined_all_strains.treefile','combined_all_strains.iqtree']:
 shutil.copy2(ROOT/'comparative_visualization/2026-11-11/sources/input_trees'/name,OUT/'data'/name)
shutil.copy2(OUT/'data/lovastatin_region_features.tsv',OUT/'data/region_features.tsv')
checks={}
checks['23_gate_rows']=len(P['gates'])==23
checks['77_unique_labels']=len({x['sample'] for x in P['genomes']})==77
checks['72_distinct_sequence_sets']=len({x['sequence_sha256'] for x in P['genomes']})==72
checks['66_core_genomes']=sum(x['core_tree'] for x in P['genomes'])==66
checks['101_lovastatin_pieces']=sum(x['segment_count'] for x in P['genomes'])==101
checks['all_18_evidence_snapshots_match_originals']=all(digest(OUT/x['snapshot'])==x['sha256']==digest(x['source']) for x in P['sources'])
figsrc=json.loads((OUT/'data/figure_source_hashes.json').read_text());checks['all_figure_sources_unchanged']=all(digest(k)==v for k,v in figsrc.items())
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(OUT/'analysis_gates.xlsx') as z:
 sheets=[n for n in z.namelist() if re.fullmatch(r'xl/worksheets/sheet\d+\.xml',n)]
 checks['xlsx_has_7_sheets']=len(sheets)==7
 checks['xlsx_has_7_filterable_tables']=len([n for n in z.namelist() if re.fullmatch(r'xl/tables/table\d+\.xml',n)])==7
 roots=[E.fromstring(z.read(n)) for n in sheets]
 checks['xlsx_no_error_cells']=not any(c.get('t')=='e' for r in roots for c in r.findall('.//s:c',ns))
 values={c.get('r'):c.find('s:v',ns).text for c in roots[0].findall('.//s:c',ns) if c.find('s:f',ns) is not None}
 counts=collections.Counter(r['quality'] for r in P['gates'])
 checks['xlsx_status_formula_caches_correct']=all(int(values[cell])==counts[status] for cell,status in [('B5','Pass'),('D5','Conditional'),('F5','Open'),('H5','Hold'),('B6','Exploratory'),('D6','Not run')])
 checks['xlsx_freeze_headers']=all(r.find('s:sheetViews/s:sheetView/s:pane',ns).get('ySplit') in ['6','8'] for r in roots)
 core_rows=roots[2].findall('.//s:row',ns)
 cells={c.get('r'):c for row in core_rows for c in row.findall('s:c',ns)}
 checks['xlsx_genome_numbers_are_numeric']=all(cells[f'{col}{row}'].get('t','n')=='n' for col in ['H','I','J','K','L','N','O'] for row in range(7,84))
pdfpages={}
for fn,want in [('01_three_gene_tree.pdf',4),('02_core_genome_tree.pdf',3),('03_lovastatin_linear.pdf',1),('04_whole_genome_circos.pdf',1),('05_marker_comparison.pdf',1)]:
 d=pymupdf.open(F/fn);pdfpages[fn]=len(d);checks[fn+'_pages']=len(d)==want
 checks[fn+'_nonempty_pages']=all(len(pg.get_text())>100 for pg in d)
 if fn=='02_core_genome_tree.pdf':
  text=d[1].get_text();checks['all_66_core_leaf_labels_present']=all(a['display_label'].replace('[TG]','[TYP-G]') in text for a in P['core_annotations'])
 if fn in ['03_lovastatin_linear.pdf','04_whole_genome_circos.pdf']:
  text=d[0].get_text();checks[fn+'_77_labels']=all(re.search(r'\|\s*'+re.escape(m['sample'])+r'(?:\s|$)',text) for m in P['genomes'])
  if fn.startswith('04'):checks['plain_circos_no_cluster_content']=not any(x in text.lower() for x in ['lovastatin','lovb','pathway'])
 d.close()
book=pymupdf.open(OUT/'stakeholder_figures.pdf');checks['figure_book_11_pages']=len(book)==11
checks['figure_book_bookmarks']=len(book.get_toc())==6
book[0].get_pixmap(matrix=pymupdf.Matrix(1,1)).save(str(Q/'book_cover.png'));book.close()
orig=pymupdf.open(ROOT/'comparative_visualization/2026-11-13/marker_provenance/combined_type_marked/kombinierter_baum_typstaemme.pdf');new=pymupdf.open(F/'01_three_gene_tree.pdf')
checks['combined_tree_existing_views_unchanged']=all(orig[i].get_text()==new[i].get_text() and orig[i].get_pixmap(matrix=pymupdf.Matrix(.12,.12)).samples==new[i].get_pixmap(matrix=pymupdf.Matrix(.12,.12)).samples for i in [0,1]);orig.close();new.close()
checks['two_main_keys_use_existing_species_palette']=json.loads((OUT/'data/sources/S08/taxon_palette.json').read_text())==P['palette']
figure_checks=json.loads((Q/'figure_validation.json').read_text())['checks'];checks.update(figure_checks)
result={'release':P['release'],'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'check_count':len(checks),'quality_status_counts':dict(counts),'pdf_page_counts':pdfpages,'resource_limits':LIMITS}
(Q/'release_validation.json').write_text(json.dumps(result,indent=2))
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
readme='''ASPERGILLUS TERREUS GROUP - STAKEHOLDER RELEASE
Version: 2026-11-17_v1

Start here
  analysis_gates.xlsx        7 sheets, 23 analysis-status / quality-gate entries
  stakeholder_figures.pdf    11 pages, bookmarked; five main figure families
  index.html                Local figure gallery with links

Individual PDFs (figures/)
  01_three_gene_tree.pdf       4 pages: latest circular + vertical combined tree,
                              full colour legend and tree-specific type-clade key
  02_core_genome_tree.pdf      3 pages: Klebsiella-style circular core tree,
                              uncompressed linear tree and metadata key
  03_lovastatin_linear.pdf     1 large page: 77 labels / 101 separate contig pieces,
                              pathway homologs and up to 20 kb flanks
  04_whole_genome_circos.pdf   1 page: all 77 labels / 72 canonical sequences,
                              assembly-wide metrics, no cluster focus
  05_marker_comparison.pdf    1 page: RF ranking / bootstrap percentiles and NMDS

PDFs retain vector graphics. SVG and PNG companions are supplied in figures/.
The 283-tip tree and 101-piece architecture figure use large pages to preserve
all labels. Zoom in or print at native size. The original two type-marked marker
pages are unchanged; the complete colour and clade keys are appended.
Species colours and the N50 scale are harmonised across the new tree/Circos set.
The Circos taxon ring is REPORTED taxon; tree species rings are REVIEWED taxon.
Their meaning remains distinct even when the named-species colour is shared.

Scientific status
  2 passed technical checks; 11 conditional; 5 open; 3 held; 1 exploratory;
  1 unrun ancestral-state workstream (five candidate approaches on its own tab).
  The core tree has no bootstrap / SH-aLRT support and remains exploratory.
  Harmonised genome/transcript BUSCO, annotation version/model integrity,
  biological-isolate alias reconciliation and functional validation remain open.
  A lovastatin homology block is not proof of intact CDS or metabolite production.
  No convincing block is not confirmed absence. No raw-read result is invented.
  Legacy secretome, SNP and functional-clustering plots are held back.
  Ancestral-state / evo-devo approaches were discussed, not computed.

Cohorts
  Core phylogeny: 66 tips, 4,661 loci, 9,012,257 alignment columns.
  Combined marker tree: 283 tips = 76 genome-derived + 207 reference entries.
  Type marks: 31 reference ex-type entries and 3 documented ex-type genomes;
  these are entry counts, not 34 independent cultures.
  Assembly displays: 77 labels, 72 distinct canonical assembly sequences.
  Alternative assembly versions are not additional independent tree isolates.

Evidence and reproducibility
  Original audit: Nov 10. Analytical corrections: Nov 10-11. Previous figure
  exports: Nov 13-14. This release redraws/packages those reviewed results;
  it does not rerun BLAST, BUSCO, MAFFT, IQ-TREE, Panaroo or phytools.
  data/package.json contains the release tables and annotations.
  data/sources/ contains 18 evidence snapshots identified in the workbook.
  data/source_manifest.json and figure_source_hashes.json retain provenance.
  data/genome_metrics.json preserves exact composition bins/contig offsets.
  data/region_features.tsv supplies exact lovastatin/flank gene identifiers.
  scripts/ retains the builders and plotting helpers. Builders use the original
  workspace release inputs and the recorded bundled runtimes, not a new inference.
  qa/release_validation.json records the technical checks and source preservation.
  Files are hashed in release_manifest.json. The companion ZIP excludes caches
  and temporary intermediate PDFs, but includes all finished outputs and QA.

Resource control / preservation
  Figure and workbook workers were limited to a separate 2 GiB Windows job and
  one CPU affinity; workbook JavaScript heap was capped at 768 MiB.
  Existing project outputs were preserved. No D: project files were modified.
  No process, input or memory limit of the Panaroo job was changed.
'''
(OUT/'README.txt').write_text(readme,encoding='utf-8')
preview={'F01':'01_three_gene_circular.png','F02':'02_core_genome_circular.png','F03':'03_lovastatin_linear.png','F04':'04_whole_genome_circos.png','F05':'05_marker_comparison.png'}
cards=[]
for ident,title,fn,coverage,purpose,caveat in P['figures']:
 cards.append(f'<article><a href="figures/{fn}"><img src="figures/{preview[ident]}" alt="{html.escape(title)}"></a><h2>{ident} · {html.escape(title)}</h2><p>{html.escape(coverage)}</p><p>{html.escape(caveat)}</p><a href="figures/{fn}">Open PDF</a></article>')
(OUT/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Terreus stakeholder release</title><style>body{max-width:1300px;margin:40px auto;padding:0 24px;font:16px Arial;color:#243d54;background:#f6f8fa}h1{font-size:30px}a{color:#087f8c}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:24px}article{background:white;padding:24px;border:1px solid #dce3e9;border-radius:10px}img{width:100%;height:360px;object-fit:contain}p{line-height:1.5}.notice{padding:18px;background:#fff3d9;margin:22px 0}</style><h1>Aspergillus terreus group</h1><p>Stakeholder release · 2026-11-17_v1</p><p><a href="analysis_gates.xlsx">Excel quality gates</a> · <a href="stakeholder_figures.pdf">Complete figure book</a> · <a href="README.txt">Scope and provenance</a></p><div class="notice">The core-genome topology remains exploratory. Lovastatin panels show homolog architecture, not proven biosynthetic function. Detailed completion gates are in the workbook.</div><main>'+''.join(cards)+'</main></html>',encoding='utf-8')
files=[p for p in OUT.rglob('*') if p.is_file() and not any(k in p.parts for k in ['node_modules','matplotlib_cache','__pycache__']) and p.name not in ['release_manifest.json','cover.pdf','marker_keys.pdf']]
manifest=[{'path':str(p.relative_to(OUT)).replace('\\','/'),'bytes':p.stat().st_size,'sha256':digest(p)} for p in files]
(OUT/'release_manifest.json').write_text(json.dumps({'release':P['release'],'files':manifest},indent=2))
zip_path=OUT.parent/'terreus_stakeholder_2026-11-17_v1.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files+[OUT/'release_manifest.json']:
  if p.name in ['cover.pdf','marker_keys.pdf']:continue
  z.write(p,OUT.name+'/'+str(p.relative_to(OUT)))
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
print('RELEASE PASS',len(checks),'checks;',zip_path,'ZIP bytes',zip_path.stat().st_size)
