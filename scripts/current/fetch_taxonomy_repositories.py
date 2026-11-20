from terreus_reference_sources import *
for url,name in [('https://api.figshare.com/v2/articles/23723277','Visagie_2024_figshare.json'),('https://zenodo.org/api/records/16605949','accepted_2025_zenodo.json')]:
    try:
        d=get(url,name,'json');print('REPOSITORY',name,d.get('title',d.get('metadata',{}).get('title')),flush=True)
        for f in d.get('files',[]):print('FILE',json.dumps(f),flush=True)
    except Exception as e:print('ERROR',name,str(e),flush=True)
get('https://pmc.ncbi.nlm.nih.gov/articles/instance/10142542/bin/pathogens-12-00579-s001.zip','pseudoalabamensis_2023_supp.zip','bytes')
