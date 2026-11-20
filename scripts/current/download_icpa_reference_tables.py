from terreus_reference_sources import *
d=json.loads((RAW/'accepted_2025_zenodo.json').read_text(encoding='utf-8'))
print('RESOLVED RECORD',d.get('id'),d.get('doi'),flush=True)
for f in d['files']:
    name='ICPA_reference_sequences_20261001.xlsx' if 'Table S1' in f['key'] else 'ICPA_accepted_species_20261001.xlsx'
    data=get(f['links']['self'],name,'bytes')
    assert hashlib.md5(data).hexdigest()==f['checksum'].split(':')[1]
    print('DOWNLOADED',name,len(data),flush=True)
save('icpa_version_lock.json',{'requested_record':'16605949','resolved_record':d.get('id'),'doi':d.get('doi'),'files':d['files']})
