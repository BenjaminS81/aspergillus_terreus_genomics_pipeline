"""Keep vector PDF masters; replace raster scales; stream true 600 dpi TIFFs."""
from pathlib import Path
import sys, math, json, hashlib, struct, zlib, re, zipfile, os, gc, time
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
LIMITS=constrain_process()
import pymupdf
import numpy as np
from PIL import Image

DPI=600;S=DPI/72;ROWS=128
T=OUT/'print_600dpi';Q=OUT/'qa'/'print_validation';A=OUT/'archive'
for p in [T,Q,A]:p.mkdir(parents=True,exist_ok=True)
pdfs=sorted(OUT.glob('*.pdf'));assert len(pdfs)==10
backup=A/'before_print_upgrade_v1.zip'
if not backup.exists():
 with zipfile.ZipFile(backup,'w',zipfile.ZIP_DEFLATED) as z:
  for p in pdfs+[OUT/'manifest.json',OUT/'README.txt']:
   if p.exists():z.write(p,p.name)

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def vector_master(path):
 before=sha(path);doc=pymupdf.open(path);page=doc[0];original_text=page.get_text()
 original_rect=tuple(page.rect);original_fonts=page.get_fonts();scales=page.get_image_info(xrefs=True)
 resource_names={im[0]:im[7] for im in page.get_images(full=True)}
 # Colour bars contain a uniform vertical gradient strip. Rebuild as vector
 # bands with the exact encoded source colours, preserving limits and labels.
 for im in scales:
  pix=pymupdf.Pixmap(doc,im['xref'])
  if pix.colorspace.n!=3:pix=pymupdf.Pixmap(pymupdf.csRGB,pix)
  rgb=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,pix.n)[:,:,:3]
  assert np.array_equal(rgb[0],rgb[-1]),'Unexpected non-gradient image'
  colors=rgb[pix.height//2];n=len(colors);ops=['q']
  for j,(r,g,b) in enumerate(colors):
   ops.append(f'{r/255:.6f} {g/255:.6f} {b/255:.6f} rg {j/n:.9f} 0 {1/n+.000001:.9f} 1 re f')
  ops.append('Q');replacement='\n'.join(ops).encode()
  name=resource_names[im['xref']].encode();pattern=rb'/'+re.escape(name)+rb'\s+Do\b';hits=0
  for xr in page.get_contents():
   content=doc.xref_stream(xr);content,count=re.subn(pattern,lambda _:replacement,content);hits+=count
   if count:doc.update_stream(xr,content)
  assert hits==1,(path.name,name,hits)
 # Remove references to unused raster resources and keep all fonts embedded.
 page.clean_contents(sanitize=True)
 meta=doc.metadata;meta.update({'title':path.stem[3:].replace('_',' ')+' | print master','subject':'All-vector single-page print master; matching RGB TIFF at 600 dpi. Print at native page size.','keywords':'vector, embedded fonts, print, 600 dpi TIFF, terreus, 2026-11-17_v2'})
 doc.set_metadata(meta)
 doc.xref_set_key(doc.pdf_catalog(),'ViewerPreferences','<< /PrintScaling /None >>')
 doc.xref_set_key(doc.pdf_catalog(),'PageLayout','/SinglePage')
 doc.xref_set_key(doc.pdf_catalog(),'OpenAction',f'[{page.xref} 0 R /XYZ null null 1]')
 tmp=path.with_suffix('.print.pending.pdf');doc.save(tmp,garbage=4,deflate=True,clean=True);doc.close()
 with pymupdf.open(tmp) as d:
  p=d[0]
  assert len(d)==1 and tuple(p.rect)==original_rect
  assert p.get_text()==original_text
  assert not p.get_image_info(),(path.name,'raster content remains')
  assert all(len(d.extract_font(f[0])[3])>0 for f in p.get_fonts())
 os.replace(tmp,path)
 return {'pdf':path.name,'previous_pdf_sha256':before,'pdf_sha256':sha(path),'native_inches':[original_rect[2]/72,original_rect[3]/72],'raster_scales_vectorised':len(scales),'all_visible_content_vector':True,'text_and_page_size_preserved':True,'fonts_embedded':True}

def write_tiff(page,path):
 """Baseline RGB Deflate TIFF, strip by strip; never allocate a full page."""
 raster_bounds=(page.rect*pymupdf.Matrix(S,S)).irect
 width=raster_bounds.width;height=raster_bounds.height
 display=page.get_displaylist();offsets=[];counts=[];strip_hashes={}
 pending=path.with_suffix('.pending.tif')
 with pending.open('w+b') as h:
  h.write(b'II'+struct.pack('<HI',42,0))
  for i,y in enumerate(range(0,height,ROWS)):
   bottom=min(y+ROWS,height)
   clip=pymupdf.Rect(0,max(0,(y-1)/S),page.rect.width,min(page.rect.height,(bottom+1)/S))
   pix=display.get_pixmap(matrix=pymupdf.Matrix(S,S),colorspace=pymupdf.csRGB,alpha=False,clip=clip)
   pixels=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,3)
   band=np.ascontiguousarray(pixels[y-pix.y:bottom-pix.y,:width,:])
   assert band.shape==(bottom-y,width,3),(i,band.shape,width,height,pix.irect)
   raw=band.tobytes();compressed=zlib.compress(raw,4);offsets.append(h.tell());counts.append(len(compressed));h.write(compressed)
   if i in [0,(height//ROWS)//2,(height-1)//ROWS]:strip_hashes[i]=hashlib.sha256(raw).hexdigest()
   del raw,compressed,pix,pixels,band
  def payload(data):
   if h.tell()%2:h.write(b'\0')
   off=h.tell();h.write(data);return off
  bits=payload(struct.pack('<HHH',8,8,8));strip_offsets=payload(struct.pack('<'+'I'*len(offsets),*offsets));strip_counts=payload(struct.pack('<'+'I'*len(counts),*counts));res=payload(struct.pack('<II',DPI,1))
  software=b'Terreus print export v2; rendered from vector PDF at 600 dpi\0';soft=payload(software)
  if h.tell()%2:h.write(b'\0')
  ifd=h.tell()
  tags=[(256,4,1,width),(257,4,1,height),(258,3,3,bits),(259,3,1,8),(262,3,1,2),(273,4,len(offsets),strip_offsets),(277,3,1,3),(278,4,1,ROWS),(279,4,len(counts),strip_counts),(282,5,1,res),(283,5,1,res),(284,3,1,1),(296,3,1,2),(305,2,len(software),soft)]
  h.write(struct.pack('<H',len(tags)))
  for tag,typ,n,value in tags:h.write(struct.pack('<HHII',tag,typ,n,value))
  h.write(struct.pack('<I',0));h.seek(4);h.write(struct.pack('<I',ifd))
 # Independently inspect the container and decode sampled strips without a full-page allocation.
 Image.MAX_IMAGE_PIXELS=None
 with Image.open(pending) as im:
  assert im.format=='TIFF' and im.size==(width,height) and im.mode=='RGB'
  assert im.info['dpi']==(600.0,600.0),im.info
  assert im.tag_v2[259]==8
 with pending.open('rb') as h:
  for i,digest in strip_hashes.items():
   h.seek(offsets[i]);raw=zlib.decompress(h.read(counts[i]));assert hashlib.sha256(raw).hexdigest()==digest
 os.replace(pending,path)
 return {'tiff':'print_600dpi/'+path.name,'dpi':DPI,'pixels':[width,height],'tiff_sha256':sha(path),'tiff_bytes':path.stat().st_size,'sampled_strips_verified':True}

records=[]
progress=Q/'progress.json'
cached={r['pdf']:r for r in json.loads(progress.read_text())} if progress.exists() else {}
for i,path in enumerate(pdfs,1):
 if path.name in cached and sha(path)==cached[path.name]['pdf_sha256'] and (OUT/cached[path.name]['tiff']).exists() and sha(OUT/cached[path.name]['tiff'])==cached[path.name]['tiff_sha256']:
  records.append(cached[path.name]);print(f'VERIFIED existing print {i}/10 {path.name}',flush=True);continue
 start=time.time();record=vector_master(path)
 with pymupdf.open(path) as d:
  page=d[0];record.update(write_tiff(page,T/(path.stem+'.tif')))
  scale=min(1,1400/page.rect.width,1600/page.rect.height)
  page.get_pixmap(matrix=pymupdf.Matrix(scale,scale)).save(str(Q/(path.stem+'_overview.png')))
  # Full-resolution excerpt for print-level inspection, as opposed to a fitted preview.
  r=page.rect
  if path.stem.startswith('01'):clip=pymupdf.Rect(r.width*.72,r.height*.42,r.width*.99,r.height*.55)
  elif path.stem.startswith('02'):clip=pymupdf.Rect(r.width*.37,r.height*.17,r.width*.70,r.height*.21)
  elif path.stem.startswith('08'):clip=pymupdf.Rect(r.width*.26,r.height*.08,r.width*.74,r.height*.105)
  elif path.stem.startswith('09'):clip=pymupdf.Rect(r.width*.025,r.height*.88,r.width*.25,r.height*.93)
  else:clip=pymupdf.Rect(r.width*.03,r.height*.08,r.width*.40,r.height*.18)
  page.get_pixmap(dpi=DPI,clip=clip).save(str(Q/(path.stem+'_600dpi_detail.png')))
 records.append(record);(Q/'progress.json').write_text(json.dumps(records,indent=2))
 print(f'PRINT {i}/10 {path.name}: {record["pixels"][0]} x {record["pixels"][1]} pixels at 600 dpi; {time.time()-start:.1f}s',flush=True);gc.collect()

# Explicit native print scaling and 100% opening zoom also apply to resumed files.
for record in records:
 path=OUT/record['pdf'];doc=pymupdf.open(path);pg=doc[0]
 doc.xref_set_key(doc.pdf_catalog(),'ViewerPreferences','<< /PrintScaling /None >>')
 doc.xref_set_key(doc.pdf_catalog(),'PageLayout','/SinglePage')
 doc.xref_set_key(doc.pdf_catalog(),'OpenAction',f'[{pg.xref} 0 R /XYZ null null 1]')
 pending=path.with_suffix('.view.pending.pdf');doc.save(pending,garbage=4,deflate=True);doc.close();os.replace(pending,path)
 record['pdf_sha256']=sha(path);record['native_print_scaling']=True;record['opening_zoom_percent']=100
manifest={'release':'2026-11-17_v2_print','policy':'single-page all-vector PDFs; lossless RGB TIFF at true 600 dpi at native page dimensions; no upsampling of preview PNGs','resource_limits':LIMITS,'files':records}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
lines=['PRINT EXPORTS - 2026-11-17_v2_print','','All ten root-level PDFs have been replaced with print masters.','Every PDF has one page, embedded fonts and fully vector visible content.','The six old raster colour scales are now vector graphics.','Trees, taxa, gene labels, data, coordinates, topology and page dimensions are preserved.','','Matching print_600dpi/*.tif files are rendered directly from these vector PDFs at','600 pixels per inch in both axes, with lossless Deflate compression and RGB colour.','These are real rendered pixels, not enlarged low-resolution PNG previews.','','IMPORTANT FOR PRINTING','Use the PDF for the sharpest scalable print; use TIFF when a raster submission is required.','Print at native size / 100%. Fitting the very large trees or stacked gene plot onto','A4/A3 makes small text physically unreadable regardless of image resolution.','PDF viewer thumbnails and the qa/ overview images are only fitted screen previews.','They do not show the resolution of the PDF or the 600-dpi print TIFF.','','Native sizes (width x height, cm) and 600-dpi dimensions:']
for r in records:lines.append(f'{r["pdf"]}: {r["native_inches"][0]*2.54:.1f} x {r["native_inches"][1]*2.54:.1f} cm; {r["pixels"][0]} x {r["pixels"][1]} px')
lines+=['','The previous PDFs are retained in archive/before_print_upgrade_v1.zip.','The Excel gate workbook remains unchanged in the preceding stakeholder release.','No sequence analyses or Panaroo processes were changed. Print workers used a','separate 2 GiB job memory cap, one CPU affinity and strip-based rendering.']
(OUT/'README.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(Q/'validation.json').write_text(json.dumps({'ten_pdfs':len(records)==10,'all_vector':all(r['all_visible_content_vector'] for r in records),'all_600_dpi':all(r['dpi']==600 for r in records),'all_strip_checks':all(r['sampled_strips_verified'] for r in records),'all_text_preserved':all(r['text_and_page_size_preserved'] for r in records)},indent=2))
for full,name in [(False,'terreus_single_page_pdfs_2026-11-17_v1.zip'),(True,'terreus_print_600dpi_2026-11-17_v2.zip')]:
 archive_path=OUT.parent/name
 with zipfile.ZipFile(archive_path,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
  for p in pdfs+[OUT/'README.txt',OUT/'manifest.json']+[OUT/'upgrade_print_exports.py']:z.write(p,p.relative_to(OUT))
  if full:
   for p in T.glob('*.tif'):z.write(p,p.relative_to(OUT))
 with zipfile.ZipFile(archive_path) as z:assert z.testzip() is None
 print('ARCHIVE',archive_path.name,archive_path.stat().st_size,flush=True)
print('PRINT EXPORT COMPLETE',flush=True)
