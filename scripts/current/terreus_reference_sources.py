from pathlib import Path
import urllib.request,urllib.parse,json,time,hashlib,datetime,xml.etree.ElementTree as ET
BASE=Path(__file__).resolve().parent
OUT=BASE/'taxonomy_reference_build'/'2026-11-10'
RAW=OUT/'sources';RAW.mkdir(parents=True,exist_ok=True)
manifest_path=OUT/'download_manifest.json'
manifest=json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.exists() else []
last=0
def get(url,name,kind='xml'):
    global last
    path=RAW/name
    if path.exists():raw=path.read_bytes()
    else:
        for attempt in range(4):
            time.sleep(max(0,.4-(time.monotonic()-last)))
            try:
                with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'TerreusTaxonomyReference/1.0'}),timeout=60) as r:raw=r.read()
                last=time.monotonic()
                if kind=='xml':ET.fromstring(raw)
                elif kind=='json':json.loads(raw)
                path.write_bytes(raw)
                manifest.append({'file':'sources/'+name,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)})
                manifest_path.write_text(json.dumps(manifest,indent=2),encoding='utf-8');break
            except Exception:
                last=time.monotonic()
                if attempt==3:raise
                time.sleep(1+attempt)
    return ET.fromstring(raw) if kind=='xml' else json.loads(raw) if kind=='json' else raw
def eu(endpoint,params,name,kind='xml'):
    return get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+endpoint+'.fcgi?'+urllib.parse.urlencode(params,doseq=True),name,kind)
def save(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
if __name__=='__main__':
    for uid in ['11003441','7426331']:
        root=eu('efetch',{'db':'pmc','id':uid,'retmode':'xml'},'PMC'+uid+'.xml')
        print('ARTICLE',uid,flush=True)
        for art in root.findall('.//article'):
            print('TITLE',''.join(art.find('./front/article-meta/title-group/article-title').itertext()),flush=True)
            for s in art.findall('.//supplementary-material'):
                print('SUPP',ET.tostring(s,encoding='unicode')[:2600],flush=True)
