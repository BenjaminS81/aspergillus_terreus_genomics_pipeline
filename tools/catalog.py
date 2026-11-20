"""Inspect source and taxon inventories without accessing external data."""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('view', choices=['summary', 'taxa', 'scripts'])
    p.add_argument('--scope', choices=['all', 'fetched', 'core', 'complete230'], default='all')
    a = p.parse_args()
    if a.view == 'summary':
        print((ROOT / 'manifests/publication_summary.json').read_text(encoding='utf-8'))
    elif a.view == 'scripts':
        rows = json.loads((ROOT / 'manifests/script_versions.json').read_text(encoding='utf-8'))
        for r in rows:
            print(r['category'], r['origin'], r['original_relative_path'], r['stored_sha256'], sep='\t')
    else:
        rows = json.loads((ROOT / 'inventory/taxon_evidence.json').read_text(encoding='utf-8'))
        field = {'fetched': 'sequence_fetched', 'core': 'used_core', 'complete230': 'used_complete230'}.get(a.scope)
        for r in rows:
            if field is None or r[field]:
                print(r['taxon_name'])

if __name__ == '__main__':
    main()
