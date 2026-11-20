from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlparse,unquote
import json,csv,zipfile
p=Path(__file__).resolve().parent/'review_reports/2026-11-10/ncbi_core_tree_metadata'
class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        self.links.extend(v for k,v in attrs if k in ['href','src'])
links=Links();links.feed((p/'ncbi_metadata_report.html').read_text(encoding='utf-8'))
missing=[]
for link in links.links:
    if not urlparse(link).scheme and not (p/unquote(link)).exists():missing.append(link)
rows=list(csv.DictReader((p/'core_tree_metadata.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
assert len(rows)==66 and not missing
assert sum(bool(r['ncbi_organism']) for r in rows)==51
assert len({r['biosample_accessions'] for r in rows if r['biosample_accessions']})==47
for fn in ['retrieval_errors.json','enrichment_errors.json','supplementary_errors.json']:assert not json.loads((p/fn).read_text(encoding='utf-8'))
with zipfile.ZipFile(p.parent/'ncbi_core_tree_metadata_bundle.zip') as z:
    assert z.testzip() is None
    assert len(z.namelist())==len([f for f in p.rglob('*') if f.is_file()])
print(json.dumps({'table_rows':66,'mapped_taxonomy_rows':51,'unique_biosamples':47,'missing_local_links':missing,'retrieval_error_logs':'empty','zip_integrity':'passed'},indent=2))
