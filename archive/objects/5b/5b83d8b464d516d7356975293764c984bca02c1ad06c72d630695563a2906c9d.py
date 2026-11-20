from pathlib import Path
import json,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis/orthology_crosscheck'
old=json.loads((P/'old_marker_genes.json').read_text());new=json.loads((P/'validated_marker_gene_mapping.json').read_text())
for d in old:
    if d['locus']=='BenA' and d.get('orthogroup')!='OG0009417':print('OLD EXCEPTION',d)
for d in new:
    if d['locus']=='BenA' and any(g.get('orthogroup')=='OG0002455' for g in d.get('genes',[])):print('NEW NEIGHBOUR',d)
for og,r in json.loads((P/'relevant_orthogroup_rows.json').read_text()).items():
    print('OG',og,'COUNTS',dict(collections.Counter(len([x for x in v.split(',') if x.strip()]) for k,v in r.items() if k!='Orthogroup')))
