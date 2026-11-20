from pathlib import Path
import json,re,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10';L=P/'reference_library';A=P/'analysis'
refs=json.loads((A/'validated_reference_metadata.json').read_text());allrefs=json.loads((L/'validated_reference_rows.json').read_text());accepted=json.loads((L/'accepted_species.json').read_text());records=json.loads((L/'historical_2011/records.json').read_text())
def norm(s):
    s=re.sub('[^A-Z0-9]','',s.upper())
    if s.startswith('DTO'):s=re.sub(r'^DTO(\d+)',lambda m:'DTO'+m[1].zfill(3),s)
    return s
bystrain={}
for r in records:
    if r['locus'] not in ['BenA','CaM']:continue
    strain=r['source_qualifiers'].get('strain',r['source_qualifiers'].get('isolate',''))
    assert strain
    key=norm(strain);d=bystrain.setdefault(key,{'reference_id':'HIST2011_'+key,'strain':strain,'deposited_organisms':[],'role':'historical_reference','type_flag':False,'accepted_type_name':'','historical_type_name':'','source_paper':'10.3114/sim.2011.69.04','section':'historical_Terrei_or_outgroup','series':'','subgenus':'','source_row':None})
    d['deposited_organisms'].append(r['organism']);assert r['locus']+'GB' not in d
    d[r['locus']+'GB']=r['accession-version'];d[r['locus']+'_analysis_sequence']=r['sequence'].upper()
for key,d in bystrain.items():
    aliases=[r for r in allrefs if norm(str(r['strain']))==key or any(str(r.get(l+'GB'))==str(d.get(l+'GB','')).split('.')[0] for l in ['BenA','CaM'] if d.get(l+'GB'))]
    types=[a for a in accepted if any(norm(x.strip())==key for x in str(a['EX_TYPUS']).split('='))]
    mapped={r['species'] for r in aliases}|{a['ID'].replace(' ','_') for a in types}
    d['species']=next(iter(mapped)) if len(mapped)==1 else sorted(set(d['deposited_organisms']))[0].replace(' ','_')
    d['analysis_species']=d['species'];d['taxonomy_status']='mapped_to_current_reference_strain' if len(mapped)==1 else 'historical_name_requires_placement_review';d['culture_alias_reference_ids']=[r['reference_id'] for r in aliases]
    for l in ['BenA','CaM','RPB2']:d.setdefault(l+'GB','');d.setdefault(l+'_analysis_sequence','')
hist=list(bystrain.values())
for l in ['BenA','CaM']:(A/(l+'_historical_2011.fasta')).write_text(''.join('>'+d['reference_id']+'\n'+d[l+'_analysis_sequence']+'\n' for d in hist if d[l+'_analysis_sequence']))
(L/'historical_2011/merged_strain_metadata.json').write_text(json.dumps(hist,indent=2));(A/'tree_reference_metadata.json').write_text(json.dumps(refs+hist,indent=2))
print('HISTORICAL STRAIN LABELS',len(hist),'MARKERS',sum(bool(d[l+'_analysis_sequence']) for d in hist for l in ['BenA','CaM']))
