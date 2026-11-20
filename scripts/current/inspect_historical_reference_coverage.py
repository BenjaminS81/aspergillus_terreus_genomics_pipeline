from pathlib import Path
import json,re,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/reference_library'
r=json.loads((P/'reference_rows.json').read_text());h=json.loads((P/'historical_2011/records.json').read_text());norm=lambda x:re.sub('[^A-Z0-9]','',str(x).upper())
for d in h:
    if d['locus']=='ITS':continue
    a=d['primary-accession'];strain=d['source_qualifiers'].get('strain',d['source_qualifiers'].get('isolate',''))
    matches=[x for x in r if x.get(d['locus']+'GB')==a or norm(x['strain'])==norm(strain)]
    print(a,d['locus'],d['organism'],strain,[(m['species'],m['strain']) for m in matches])
