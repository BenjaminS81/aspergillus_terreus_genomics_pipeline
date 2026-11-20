from pathlib import Path
import json,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10'
f=P/'analysis/genomic_markers/locus_calls.json'
if f.exists():
    calls=json.loads(f.read_text());print('CALLS',collections.Counter((d['locus'],d['status']) for d in calls))
    for d in calls:
        if d['status']!='single_strong_candidate':print('FLAG',d['sample'],d['locus'],d['status'],d['candidates'])
    for l in ['BenA','CaM','RPB2']:
        x=[d for d in calls if d['locus']==l and d.get('sequence')];print('RANGE',l,len(x),min(d['weighted_identity_percent'] for d in x),min(d['query_coverage_percent'] for d in x))
for d in json.loads((P/'reference_library/accepted_species.json').read_text()):
    if 'flavipes' in d['ID']:print('OUTGROUP',d['ID'],d['EX_TYPUS'])
