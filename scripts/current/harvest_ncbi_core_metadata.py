"""Retrieve public NCBI metadata for exactly the previously audited core-tree tips.
Writes only into the C: review workspace; source data are read-only.
"""
from pathlib import Path
import json, re, csv, io, time, hashlib, datetime, urllib.request, urllib.parse
import xml.etree.ElementTree as ET
import argparse

BASE=Path(__file__).resolve().parent
ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
OUT=BASE/'review_reports'/'2026-11-10'/'ncbi_core_tree_metadata'
RAW=OUT/'raw'; RAW.mkdir(parents=True,exist_ok=True)
E=json.loads((BASE/'review_reports/2026-11-10/audit_evidence.json').read_text())
MANIFEST=OUT/'download_manifest.json'
downloads=json.loads(MANIFEST.read_text()) if MANIFEST.exists() else []
last=0

def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
def get(url,name,kind='json'):
    global last
    p=RAW/name
    if p.exists():data=p.read_bytes()
    else:
        for attempt in range(4):
            time.sleep(max(0,.38-(time.monotonic()-last)))
            try:
                req=urllib.request.Request(url,headers={'User-Agent':'TerreusMetadataReview/1.0 (public research metadata retrieval)'})
                with urllib.request.urlopen(req,timeout=50) as r:data=r.read()
                last=time.monotonic()
                if kind=='json':json.loads(data)
                if kind=='xml':ET.fromstring(data)
                p.write_bytes(data)
                downloads.append({'file':'raw/'+name,'url':url,'retrieved_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
                MANIFEST.write_text(json.dumps(downloads,indent=2),encoding='utf-8')
                break
            except Exception as exc:
                last=time.monotonic()
                if attempt==3:raise
                time.sleep(1+attempt)
    if kind=='json':return json.loads(data)
    if kind=='xml':return ET.fromstring(data)
    return data.decode('utf-8-sig',errors='replace')
def eu(endpoint,params,name,kind='json'):
    return get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+endpoint+'.fcgi?'+urllib.parse.urlencode(params,doseq=True),name,kind)
def search(db,terms,label):
    q=' OR '.join('"'+t+'"[All Fields]' for t in sorted(set(terms)))
    if not terms:return []
    result=eu('esearch',{'db':db,'term':q,'retmode':'json','retmax':10000},label+'_search.json')
    assert int(result['esearchresult']['count'])<=10000
    return result['esearchresult']['idlist']
def text(el,path):
    n=el.find(path);return ''.join(n.itertext()).strip() if n is not None else ''
def xmlstr(el):return ET.tostring(el,encoding='unicode')
def ids_in(obj,pattern):return sorted(set(re.findall(pattern,json.dumps(obj))))

def discover():
    rows=[]
    for tip,name in zip(E['trees']['genome']['tips'],E['trees']['genome']['normalized_tips']):
        r={'tree_tip':tip,'sample_id':name,'tree_path':E['trees']['genome']['path'],'tree_sha256':E['trees']['genome']['sha256'],
           'run_accessions':re.findall(r'[SED]RR\d+',name),'assembly_accessions':[],'local_sources':[],'local_genbank_references':[]}
        folders=list((ROOT/'data/assemblies').glob('*/'+name))
        for folder in folders:
            for fa in folder.glob('GC[AF]_*_genomic.fna'):
                r['assembly_accessions']+=re.findall(r'GC[AF]_\d+\.\d+',fa.name)
                r['local_sources'].append(str(fa))
            gb=folder/'genomic.gbff'
            if gb.exists():
                lines=[]
                with gb.open(errors='replace') as f:
                    for line in f:
                        if line.startswith('ORIGIN'):break
                        lines.append(line)
                head=''.join(lines)
                r['local_genbank_header']=head
                r['local_genbank_pubmed']=re.findall(r'PUBMED\s+(\d+)',head)
                r['local_bioprojects']=sorted(set(re.findall(r'PRJ[END][AB]\d+',head)))
                r['local_biosamples']=sorted(set(re.findall(r'SAM[NED][A-Z]?\d+',head)))
                r['local_sources'].append(str(gb))
        if r['run_accessions']:r['mapping_status']='public_run_id_in_tree_tip'
        elif r['assembly_accessions']:r['mapping_status']='public_assembly_id_in_matching_local_folder'
        else:r['mapping_status']='local_isolate_no_public_accession_established'
        r['mapping_caveat']='Accession association is based on tree label and local source naming; exact historical annotation-to-genome identity remains a separate release gate.'
        rows.append(r)
    save('tip_source_mapping.json',rows)
    print('DISCOVER',collections_summary(rows),flush=True)
    return rows
def collections_summary(rows):
    return {'tips':len(rows),'with_run':sum(bool(r['run_accessions']) for r in rows),'with_assembly':sum(bool(r['assembly_accessions']) for r in rows),'local_only':sum(not(r['run_accessions'] or r['assembly_accessions']) for r in rows)}

def collect():
    rows=discover(); errors=[]; assemblies=[]
    for acc in sorted({a for r in rows for a in r['assembly_accessions']}):
        try:
            d=get('https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+acc+'/dataset_report',acc+'_dataset_report.json')
            assemblies+=d.get('reports',[])
            print('ASSEMBLY',acc,flush=True)
        except Exception as exc:errors.append({'stage':'assembly','accession':acc,'error':str(exc)})
    save('assembly_reports.json',assemblies)
    runs=sorted({a for r in rows for a in r['run_accessions']})
    sraids=search('sra',runs,'sra_tree_runs')
    sra=eu('efetch',{'db':'sra','id':','.join(sraids),'retmode':'xml'},'sra_tree_runs.xml','xml')
    runinfo=eu('efetch',{'db':'sra','id':','.join(sraids),'rettype':'runinfo','retmode':'text'},'sra_tree_runs_runinfo.csv','text')
    runrows=list(csv.DictReader(io.StringIO(runinfo)));save('runinfo_records.json',runrows)
    packages=[]
    for pkg in sra.findall('.//EXPERIMENT_PACKAGE'):
        packages.append({'run_accessions':[r.get('accession') for r in pkg.findall('.//RUN')],
                         'experiment_accession':pkg.find('EXPERIMENT').get('accession'),
                         'biosamples':ids_in(xmlstr(pkg),r'SAM[NED][A-Z]?\d+'),
                         'bioprojects':ids_in(xmlstr(pkg),r'PRJ[END][AB]\d+'),
                         'study_title':text(pkg,'.//STUDY_TITLE'),'study_abstract':text(pkg,'.//STUDY_ABSTRACT'),
                         'sample_title':text(pkg,'.//SAMPLE/TITLE'),
                         'sample_attributes':[{'name':text(a,'TAG'),'value':text(a,'VALUE')} for a in pkg.findall('.//SAMPLE_ATTRIBUTE')],
                         'library_strategy':text(pkg,'.//LIBRARY_STRATEGY'),'library_source':text(pkg,'.//LIBRARY_SOURCE'),
                         'library_selection':text(pkg,'.//LIBRARY_SELECTION'),
                         'study_links':[{'db':text(a,'DB'),'id':text(a,'ID')} for a in pkg.findall('.//XREF_LINK')]})
    save('sra_packages_parsed.json',packages)
    print('SRA',len(sraids),'packages',len(packages),'runinfo rows',len(runrows),flush=True)
    bios=ids_in([assemblies,packages,rows],r'SAM[NED][A-Z]?\d+')
    projects=ids_in([assemblies,packages,rows],r'PRJ[END][AB]\d+')
    indices={}
    for db,accessions in [('biosample',bios),('bioproject',projects),('assembly',sorted({a for r in rows for a in r['assembly_accessions']}))]:
        uids=search(db,accessions,db+'_all');indices[db]={'requested_accessions':accessions,'uids':uids}
        if db=='assembly':
            eu('esummary',{'db':db,'id':','.join(uids),'retmode':'json'},db+'_all_summary.json')
        else:
            result=eu('efetch',{'db':db,'id':','.join(uids),'retmode':'xml'},db+'_all.xml','xml')
        print(db,len(accessions),'accessions',len(uids),'records',flush=True)
    indices['sra']={'requested_accessions':runs,'uids':sraids}
    save('entrez_indices.json',indices)
    save('retrieval_errors.json',errors)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--discover',action='store_true');args=ap.parse_args()
    discover() if args.discover else collect()
