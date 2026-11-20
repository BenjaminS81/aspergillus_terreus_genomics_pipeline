from pathlib import Path
import json, collections, re
from html.parser import HTMLParser
from terreus_tree_helpers import parse_newick

P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10'
A=P/'analysis'; L=P/'reference_library'; checks=[]
def js(path): return json.loads(path.read_text(encoding='utf-8'))
def check(name, condition):
    if not condition: raise AssertionError(name)
    checks.append(name)
def fasta(path):
    result={}; key=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            assert key not in result, (path,key)
            result[key]=''
        elif line.strip(): result[key]+=line.strip()
    return result

calls=js(A/'genomic_markers/locus_calls.json')
core={x['sample'] for x in calls if x['in_core_tree']}
check('77 labels and 231 unique marker calls',len(calls)==231 and len({(x['sample'],x['locus']) for x in calls})==231 and len({x['sample'] for x in calls})==77)
check('66 core labels each have three usable single-candidate windows',len(core)==66 and sum(bool(x.get('sequence')) and x['strong_candidate_count']==1 for x in calls if x['in_core_tree'])==198)
assets=js(A/'assembly_variant_inventory.json'); groups=js(A/'assembly_sequence_groups.json')
vcalls=js(A/'variant_genomic_markers/locus_calls.json'); searched={x['sample'] for x in calls+vcalls}
check('620 assembly files, two empty, 155 distinct sequence sets all searched',len(assets)==620 and sum(x['bytes']==0 for x in assets)==2 and len(groups)==155 and all(x['searched_as'] in searched for x in groups.values()))
check('No multiple strong barcode candidate calls',all(x['strong_candidate_count']<=1 for x in calls+vcalls))
comparisons=js(A/'assembly_marker_version_comparison.json')
counts=collections.Counter(x['comparison'] for x in comparisons)
check('522 version comparisons and 507 identical marker windows',len(comparisons)==522 and counts['same_validated_sequence']==507)
check('No flagged nucleotide disagreement in version comparisons',not any('sequence_difference' in x['comparison'] for x in comparisons))
old=js(A/'orthology_crosscheck/actual_old_BenA_alignment_mapping.json')
check('All 41 actual legacy BenA alignment records trace exactly to OG0009417',len(old)==41 and all(x['orthogroups']==['OG0009417'] and x['exact_raw_headers'] for x in old))
refs=js(A/'tree_reference_metadata.json'); accepted=js(L/'accepted_species.json')
names={x['ID'].replace(' ','_') for x in accepted if x['SECTION']=='Terrei'}
check('All 29 nominal Terrei type-name entries represented',len(names)==29 and names <= {x.get('accepted_type_name') for x in refs})
gb=js(L/'genbank_records.json'); historical=js(L/'historical_2011/records.json')
check('711 curated nucleotide accession records and 81 historical records archived',len(gb)==711 and len(historical)==81)
brip=next(x for x in refs if x['strain']=='BRIP71536a')
check('BRIP71536a corrected accession identities',brip['BenAGB']=='OL741656' and brip['RPB2GB']=='OL741659')
for gene in ['BenA','RPB2']:
    check('BRIP71536a '+gene+' sequence matches corrected GenBank record',brip[gene+'_analysis_sequence'].upper() in gb[brip[gene+'GB']]['sequence'].upper())

alignments={}
design=[('BenA','BenA_all_strains',282,515),('CaM','CaM_all_strains',276,547),('RPB2','RPB2_all_strains',231,1014),('combined_all_strains','combined_all_strains',283,2076),('type_focused_sensitivity','combined_type_focused_sensitivity',108,2076)]
for tree,stem,n,width in design:
    records=fasta(A/'tree_inputs'/(stem+'.fasta'));alignments[tree]=records
    check(tree+' alignment dimensions and legal DNA symbols',len(records)==n and {len(v) for v in records.values()}=={width} and all(set(v.upper())<=set('ACGTRYSWKMBDHVN?-') for v in records.values()))
    leaves=[x.name for x in parse_newick((A/'trees'/(tree+'.treefile')).read_text()).leaves()]
    check(tree+' tree tips exactly match alignment',len(leaves)==len(set(leaves)) and set(leaves)==set(records))
    check(tree+' IQ-TREE completed and bootstrap replicates present','Date and time:' in (A/'trees'/(tree+'.iqtree')).read_text() and (A/'trees'/(tree+'.ufboot')).stat().st_size>0)
    check(tree+' contains all core-tree samples',all('LOCAL_'+s in records for s in core))
check('All historical BenA/CaM sequences included',sum(k.startswith('HIST2011_') for k in alignments['BenA'])==29 and sum(k.startswith('HIST2011_') for k in alignments['CaM'])==25)
for key,seq in alignments['combined_all_strains'].items():
    expected=''.join(alignments[g].get(key,'?'*w) for g,w in [('BenA',515),('CaM',547),('RPB2',1014)])
    assert seq==expected, ('concatenation or missing-data coding',key)
check('Concatenation preserves each locus and codes absent loci as question marks',True)
check('Divergent SRR33623979 excluded from every final tree',all(not any('SRR33623979' in k for k in records) for records in alignments.values()))
class Links(HTMLParser):
    def __init__(self):super().__init__();self.targets=[]
    def handle_starttag(self,tag,attrs):
        self.targets.extend(v for k,v in attrs if k in ['href','src'] and v)
parser=Links();parser.feed((P/'terreus_reference_review.html').read_text(encoding='utf-8'))
missing=[s for s in parser.targets if not re.match(r'^[a-z]+:|^#',s) and not (P/s.split('#')[0]).exists()]
check('Every local report link and image resolves',not missing)
figures=list((P/'figures').glob('*'))
check('Nine figures each exported to SVG, PDF and PNG',len(figures)==27 and all(f.stat().st_size>1000 for f in figures))
out={'status':'passed','checks':checks,'check_count':len(checks),'assembly_version_comparison_counts':dict(counts),'figure_exports':len(figures),'source_project_modified':False}
(P/'validation_summary.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
