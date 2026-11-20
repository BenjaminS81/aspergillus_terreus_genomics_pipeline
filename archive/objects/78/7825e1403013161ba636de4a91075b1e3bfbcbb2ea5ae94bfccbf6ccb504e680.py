from pathlib import Path
import json,collections,csv
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
groups=json.loads((P/'assembly_sequence_groups.json').read_text());canonical=json.loads((P/'genomic_markers/locus_calls.json').read_text());variants=json.loads((P/'variant_genomic_markers/locus_calls.json').read_text());lookup={(d['sample'],d['locus']):d for d in canonical+variants};records=[]
for sha,group in groups.items():
    for sample in group['samples']:
        for locus in ['BenA','CaM','RPB2']:
            c=lookup.get((sample,locus));v=lookup[group['searched_as'],locus];a=(c or {}).get('sequence','');b=v.get('sequence','')
            status='same_validated_sequence' if a and b and a==b else 'compatible_terminal_window_difference' if a and b and (a in b or b in a) else 'sequence_difference_review' if a and b else 'no_validated_sequence_for_comparison'
            records.append({'sample':sample,'locus':locus,'assembly_sequence_sha256':sha,'searched_as':group['searched_as'],'comparison':status,'canonical_status':(c or {}).get('status',''),'variant_status':v['status'],'variant_strong_candidates':v['strong_candidate_count'],'canonical_bases':len(a),'variant_bases':len(b),'source_paths':';'.join(group['files'])})
(P/'assembly_marker_version_comparison.json').write_text(json.dumps(records,indent=2))
with (P/'assembly_marker_version_comparison.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
print('VARIANT COMPARISON',collections.Counter(d['comparison'] for d in records))
print('COPY CANDIDATES',collections.Counter(d['strong_candidate_count'] for d in canonical+variants))
