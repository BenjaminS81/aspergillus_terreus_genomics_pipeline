from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlparse, unquote
import json

folder=Path(__file__).resolve().parent/'review_reports'/'2026-11-10'
class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[];self.ids=[]
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ('href','src'):self.links.append(v)
            if k=='id':self.ids.append(v)
p=Links();p.feed((folder/'terreus_deep_review.html').read_text(encoding='utf-8'))
missing=[];checked=0
for link in p.links:
    u=urlparse(link)
    if u.scheme in ('http','https'):continue
    if link.startswith('#'):
        if link[1:] not in p.ids:missing.append(link)
        continue
    target=Path(unquote(u.path).lstrip('/')) if u.scheme=='file' else folder/unquote(u.path)
    checked+=1
    if not target.exists():missing.append(str(target))
e=json.loads((folder/'audit_evidence.json').read_text())
g=json.loads((folder/'global_status_and_gates.json').read_text())
assert sum(int(v) for v in e['orthogroups']['occupancy_histogram'].values())==14964
assert len(g['stages'])==19 and len(g['gates'])==6
print(json.dumps({'local_links_checked':checked,'missing_links':missing,'report_bytes':(folder/'terreus_deep_review.html').stat().st_size,'stages':len(g['stages'])},indent=2))
assert not missing
