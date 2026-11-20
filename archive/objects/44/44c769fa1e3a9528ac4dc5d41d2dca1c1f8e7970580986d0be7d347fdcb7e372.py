from pathlib import Path
import json
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/reference_library'
r={d['reference_id']:d for d in json.loads((P/'reference_rows.json').read_text())};g=json.loads((P/'genbank_records.json').read_text());v=json.loads((P/'genbank_reconciliation.json').read_text())
for d in v:
    if d.get('resolved_status')=='substantial_curated_genbank_difference':
        seq=r[d['reference_id']][d['locus']+'_sequence'];record=g[d['accession']]
        print(d['reference_id'],d['locus'],d['accession'],record['definition'],record['length'])
        print('UNMATCHED',[(i,seq[i-1]) for i in d['unmatched_query_positions'][:30]])
        print('STARTS',seq[:55],record['sequence'][:55])
