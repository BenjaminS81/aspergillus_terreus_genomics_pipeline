from pathlib import Path
import json,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10'
L=P/'reference_library';rows=json.loads((L/'reference_rows.json').read_text());gb=json.loads((L/'genbank_records.json').read_text());rec={(d['reference_id'],d['locus']):d for d in json.loads((L/'genbank_reconciliation.json').read_text())}
changes=[]
for d in rows:
    for l in ['BenA','CaM','RPB2']:
        seq=d[l+'_sequence'];v=rec[d['reference_id'],l];status=v.get('resolved_status',v['status']);d[l+'_original_accession']=d[l+'GB']
        if d['source_row']==3494 and l in ['BenA','RPB2']:
            corrected='OL741656' if l=='BenA' else 'OL741659';full=gb[corrected]['sequence'].upper()
            assert seq in full,(l,corrected)
            changes.append({'reference_id':d['reference_id'],'locus':l,'action':'correct_accession_gene_swap','original_accession':d[l+'GB'],'corrected_accession':corrected,'evidence':'curated sequence exactly contained in alternate accession; GenBank gene annotation agrees'})
            d[l+'GB']=corrected;status='exact_subsequence_after_accession_swap_correction'
        elif status in ['minor_curated_genbank_difference','substantial_curated_genbank_difference']:
            positions=set(v['unmatched_query_positions'])|{m['query_position'] for m in v['mismatches']}
            seq=''.join('N' if i+1 in positions else c for i,c in enumerate(seq))
            changes.append({'reference_id':d['reference_id'],'locus':l,'action':'mask_curated_bases_not_supported_by_current_genbank','positions_1based':sorted(positions),'masked_bases':len(positions),'accession':d[l+'GB']})
            status='conservative_mask_of_discrepancies'
        d[l+'_analysis_sequence']=seq;d[l+'_verification_status']=status
(L/'validated_reference_rows.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False));(L/'reference_corrections.json').write_text(json.dumps(changes,indent=2))
selected=json.loads((P/'analysis/reference_metadata.json').read_text());byid={d['reference_id']:d for d in rows}
for d in selected:
    orig=byid[d['reference_id']]
    for l in ['BenA','CaM','RPB2']:
        for suffix in ['GB','_original_accession','_analysis_sequence','_verification_status']:d[l+suffix]=orig[l+suffix]
    d['analysis_species']=d['species']
    if d['source_row']==3494:d['historical_type_name']='Aspergillus_marneyi'
    elif d['accepted_type_name'] and d['accepted_type_name']!=d['species']:d['historical_type_name']=d['accepted_type_name']
    else:d['historical_type_name']=''
(P/'analysis/validated_reference_metadata.json').write_text(json.dumps(selected,indent=2,ensure_ascii=False))
print('Validated references',len(rows),'CORRECTIONS',len(changes))
