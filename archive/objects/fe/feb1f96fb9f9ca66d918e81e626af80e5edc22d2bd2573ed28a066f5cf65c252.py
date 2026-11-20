from pathlib import Path
import sys,json,subprocess,gc,hashlib
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
limits=constrain_process()
import fitz
from PIL import Image
OUT=ROOT/'stakeholder_release/2026-11-20_v1_type_monophyly'
poppler=Path(r'C:\Users\q287js\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
stems=json.loads((OUT/'data/figure_stems.json').read_text());checks=[]
for stem in stems:
 pdf=OUT/'figures'/(stem+'.pdf');doc=fitz.open(pdf);assert len(doc)==1
 page=doc[0];outside=[]
 for b in page.get_text('dict')['blocks']:
  for line in b.get('lines',[]):
   for s in line['spans']:
    x0,y0,x1,y1=s['bbox']
    if x0<0 or y0<0 or x1>page.rect.width+.1 or y1>page.rect.height+.1:outside.append(s['text'])
 assert not outside,(stem,outside)
 subprocess.run([str(poppler),'-r','55','-singlefile','-png',str(pdf),str(OUT/'qa'/stem)],check=True,capture_output=True)
 if '--preview-only' not in sys.argv:
  pix=page.get_pixmap(dpi=300,alpha=False);pix.save(str(OUT/'figures'/(stem+'.png')))
  width,height=pix.width,pix.height;del pix;gc.collect()
 else:width=height=None
 checks.append({'figure':stem,'pages':1,'text_outside_page':outside,'width_points':page.rect.width,'height_points':page.rect.height,'png_dpi':300 if width else None,'png_pixels':[width,height],'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()})
 if stem in ['three_gene230_type_monophyly_linear','core66_type_monophyly_linear']:
  for part,(a,b) in enumerate([(0,.34),(.33,.67),(.66,1)]):
   clip=fitz.Rect(page.rect.width*.36,page.rect.height*a,page.rect.width,page.rect.height*b)
   pix=page.get_pixmap(dpi=100,alpha=False,clip=clip);pix.save(str(OUT/'qa'/(stem+'_detail_'+str(part)+'.png')));del pix
 doc.close();gc.collect();print('Verified/rendered',stem,flush=True)
(OUT/'qa/automated_checks.json').write_text(json.dumps({'resource_limits':limits,'checks':checks},indent=2))
manifest=[{'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted((OUT/'figures').iterdir())]
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('QA DONE')
