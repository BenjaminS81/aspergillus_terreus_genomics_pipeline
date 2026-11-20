from pathlib import Path
import zipfile,xml.etree.ElementTree as ET
from pypdf import PdfReader
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/sources'
for name in ['Visagie_2024_review','Brazil_caves_2024']:
    reader=PdfReader(P/(name+'.pdf'));texts=[p.extract_text(extraction_mode='layout') for p in reader.pages]
    (P/(name+'.txt')).write_text('\n'.join('PAGE '+str(i+1)+'\n'+t for i,t in enumerate(texts)),encoding='utf-8')
    print(name,len(texts),'pages; Terrei pages',[i+1 for i,t in enumerate(texts) if 'Terrei' in t])
with zipfile.ZipFile(P/'Visagie_2024_supplements.zip') as z:
    print('SUPPLEMENT FILES',[(i.filename,i.file_size) for i in z.infolist()])
    dest=P/'Visagie_2024_supplements';dest.mkdir(exist_ok=True)
    for i in z.infolist():
        if i.is_dir():continue
        target=(dest/i.filename).resolve();assert dest.resolve() in target.parents
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(i))
r=ET.parse(P/'PMC12786731.xml').getroot()
for a in r.findall('.//article'):
    print('2025 TITLE',''.join(a.find('front/article-meta/title-group/article-title').itertext()),'bytes',len(ET.tostring(a)))
    for x in a.findall('.//p'):
        t=''.join(x.itertext())
        if 'sect. Terrei' in t:print(t[:1800])
