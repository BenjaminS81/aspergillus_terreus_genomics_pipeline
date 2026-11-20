from pathlib import Path
from difflib import SequenceMatcher
import json,collections,csv
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/reference_library'
rows=json.loads((P/'reference_rows.json').read_text());gb=json.loads((P/'genbank_records.json').read_text());val=json.loads((P/'genbank_validation.json').read_text());byid={d['reference_id']:d for d in rows}
IUP={k:set(v) for k,v in {'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','K':'GT','M':'AC','S':'CG','W':'AT','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}.items()}
for v in val:
    if v['status']!='sequence_discrepancy_review':continue
    d=byid[v['reference_id']];l=v['locus'];q=d[l+'_sequence'];s=gb[v['accession'].split('.')[0]]['sequence'].upper()
    rev=s.translate(str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN'))[::-1]
    orientations=[(sum(b.size for b in SequenceMatcher(None,x,q,autojunk=False).get_matching_blocks()),label,x) for label,x in [('forward',s),('reverse_complement',rev)]]
    _,orientation,s=max(orientations,key=lambda t:t[0]);v['genbank_alignment_orientation']=orientation
    blocks=SequenceMatcher(None,s,q,autojunk=False).get_opcodes();indices=[]
    for tag,a,b,c,e in blocks:
        if tag in ['equal','replace'] and b-a==e-c:indices.extend(zip(range(a,b),range(c,e)))
        elif tag!='delete':indices.extend((-1,i) for i in range(c,e))
    overlap=0;compatible=0;mismatch=[];unmatched=[]
    for si,qi in indices:
        if qi<0:continue
        if si<0:unmatched.append(int(qi)+1);continue
        overlap+=1
        if IUP[s[si]]&IUP[q[qi]]:compatible+=1
        else:mismatch.append({'query_position':int(qi)+1,'genbank_position':int(si)+1,'curated':q[qi],'genbank':s[si]})
    v.update({'overlap_bases':overlap,'iupac_compatible_bases':compatible,'query_length':len(q),'mismatches':mismatch,'unmatched_query_positions':unmatched,'compatible_percent_query':round(100*compatible/len(q),3)})
    v['resolved_status']='iupac_compatible_subsequence' if compatible==len(q) else 'minor_curated_genbank_difference' if compatible/len(q)>=.98 else 'substantial_curated_genbank_difference'
    print(v['reference_id'],l,v['resolved_status'],v['compatible_percent_query'],'mismatches',len(mismatch),'unmatched',len(unmatched))
for v in val:
    if v['accession'].startswith(('GCA','GCF')) and v['sequence_length']:v['resolved_status']='assembly_derived_curated_sequence_not_nuccore_accession'
(P/'genbank_reconciliation.json').write_text(json.dumps(val,indent=2))
print('COUNTS',collections.Counter(v.get('resolved_status',v['status']) for v in val))
