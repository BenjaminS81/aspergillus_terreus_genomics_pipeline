from terreus_reference_sources import *
d=get('https://rest.uniprot.org/uniprotkb/Q2U2U3.json','UniProt_Q2U2U3.json','json')
print(json.dumps({k:d.get(k) for k in ['primaryAccession','uniProtkbId','organism','genes','proteinDescription']},indent=2))
