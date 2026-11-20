from pathlib import Path
import json
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
d=json.loads((P/'nearest_type_concordance.json').read_text())
for x in d:
    if x['id'] in ['LOCAL_26','LOCAL_11','LOCAL_48','LOCAL_IFO6365','LOCAL_TN-484','LOCAL_SRR25636141','LOCAL_R22','LOCAL_NIH2624','LOCAL_AT2022']:
        print(x['id'],[(l,y['best_reference_species'],y['best_identity_percent']) for l,y in x['per_locus'].items()])
