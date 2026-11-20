from pathlib import Path
import xml.etree.ElementTree as ET,re,json
P=Path('taxonomy_reference_build/2026-11-10/sources')
for uid in ['7426331','3161753','10142542']:
    r=ET.parse(P/f'PMC{uid}.xml').getroot();print('\nARTICLE',uid)
    for tab in r.findall('.//table-wrap'):
        title=' '.join(''.join(x.itertext()) for x in tab.findall('label')+tab.findall('caption'))
        print('TABLE',tab.get('id'),title[:200],len(tab.findall('.//tr')))
        if uid!='7426331':
            for tr in tab.findall('.//tr')[:8]:print(' | '.join(''.join(c.itertext()).strip() for c in tr))
    for cap in r.findall('.//fig/caption'):
        t=''.join(cap.itertext())
        if 'Terrei' in t or uid!='7426331':print('FIG',t[:1800])
    for p in r.findall('.//p'):
        t=''.join(p.itertext())
        if ('outgroup' in t.lower() or 'treebase' in t.lower()) and (uid!='7426331' or 'Terrei' in t):print('METHOD',t[:1500])
