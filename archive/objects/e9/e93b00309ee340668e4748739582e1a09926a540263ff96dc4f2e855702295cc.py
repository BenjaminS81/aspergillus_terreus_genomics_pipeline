"""Stream native-size 600 dpi print images without full-page RAM allocation."""
from pathlib import Path
import sys, ast, json, math, hashlib, struct, zlib, os, zipfile, time
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
import pymupdf
import numpy as np
from PIL import Image
DPI=600; S=DPI/72; ROWS=128
# Freeze only the two reviewed utility definitions; never execute the older release script.
source=ROOT/'stakeholder_release/2026-11-17_v1_single_page_pdfs/upgrade_print_exports.py'
tree=ast.parse(source.read_text(encoding='utf-8'))
selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in {'sha','write_tiff'}]
assert len(selected)==2
exec(compile(ast.Module(body=selected,type_ignores=[]),str(source),'exec'))
(OUT/'scripts/tiff_writer_frozen.py').write_text('\n\n'.join(ast.unparse(n) for n in selected)+'\n',encoding='utf-8')
records=[]
for path in sorted(OUT.glob('*.pdf')):
 start=time.time()
 with pymupdf.open(path) as doc:
  assert len(doc)==1
  page=doc[0]
  assert not page.get_image_info()
  assert all(doc.extract_font(f[0])[3] for f in page.get_fonts())
  record={'pdf':path.name,'pdf_sha256':sha(path),'all_vector':True,'fonts_embedded':True,'native_inches':[page.rect.width/72,page.rect.height/72]}
  record.update(write_tiff(page,OUT/'print_600dpi'/(path.stem+'.tif')))
 records.append(record)
 print(path.name,record['pixels'],f'{time.time()-start:.1f}s',flush=True)
assert len(records)==4
(OUT/'manifest.json').write_text(json.dumps({'release':'2026-11-17_v3_taxonomy_review','exported':'2026-11-18','resource_limits':LIMITS,'files':records},indent=2))
(OUT/'README.txt').write_text('''TERREUS TAXONOMY FIGURE REVIEW - v3

Four single-page, all-vector PDFs with embedded fonts, matching SVGs and
losslessly compressed RGB TIFF images rendered directly at 600 dpi.
The qa/ PNGs are screen previews, not print exports.

01 Combined BenA + CaM + RPB2 tree, vertical phylogram.
02 The same tree, circular cladogram (topology rather than branch-length scale).
03 Matched 66-genome RF/NMDS comparison; unchanged distances and coordinates.
04 ONE consolidated legend: all colours, symbols, tracks, type-clade codes
   and statistical definitions used across these figures.

Both trees retain all 283 tips: 76 genome-derived and 207 reference entries.
31 reference ex-type entries and 3 documented genome type anchors are marked.
Counts are entries, not independent biological cultures. Existing unresolved
names and withheld conflicting type claims are preserved.
All annotation-track colours are shared by both tree renderings and audited
against the legend. English labels and aligned track/legend headings replace
the previous layout. Editorial interpretation has been removed from RF/NMDS.
BenA + CaM + RPB2 is the taxonomic reference framework for this dataset.

Print the PDF at native size for readable tip labels; fit-to-A4 reduces text.
The TIFFs are true 600 dpi at the PDF native dimensions, not enlarged previews.
Sources, topology, branch lengths, evidence and numerical results are unchanged.
No upstream scientific computation or Panaroo process was modified.
Export workers are separately capped at 2 GiB and one CPU.

Provenance and validation: data/source_hashes.json, data/tip_colour_audit.tsv,
qa/figure_validation.json, qa/legend_text_overlap_check.json and manifest.json.
''',encoding='utf-8')
archive=OUT.parent/'terreus_taxonomy_review_2026-11-17_v3.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
 for p in OUT.rglob('*'):
  if p.is_file() and 'matplotlib_cache' not in p.parts and p.suffix not in {'.png','.pyc'}:
   z.write(p,p.relative_to(OUT))
with zipfile.ZipFile(archive) as z: assert z.testzip() is None
print('COMPLETE',archive,archive.stat().st_size,flush=True)
