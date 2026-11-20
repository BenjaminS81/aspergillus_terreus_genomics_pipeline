from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,datetime
P=Path(__file__).resolve().parent/'comparative_visualization/2026-11-11/sources';P.mkdir(parents=True,exist_ok=True)
query='organism_id:341663 AND (gene:ATEG_09960 OR gene:ATEG_09961 OR gene:ATEG_09962 OR gene:ATEG_09963 OR gene:ATEG_09964 OR gene:ATEG_09965 OR gene:ATEG_09966 OR gene:ATEG_09967 OR gene:ATEG_09968)'
url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':query,'format':'json','size':50})
urls={'uniprot_lov_cluster_nine.json':url,'MIBiG_BGC0000088.html':'https://mibig.secondarymetabolites.org/repository/BGC0000088.4/index.html'}
manifest=[]
for name,url in urls.items():
    f=P/name
    if not f.exists():f.write_bytes(urllib.request.urlopen(url,timeout=40).read())
    manifest.append({'file':name,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(P/'download_manifest.json').write_text(json.dumps(manifest,indent=2))
for d in json.loads((P/'uniprot_lov_cluster_nine.json').read_text())['results']:print(d['primaryAccession'],d.get('genes'),d['sequence']['length'])
