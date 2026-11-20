from pathlib import Path
import json
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';raw=json.loads((BASE/'review_reports/2026-11-10/ncbi_core_tree_metadata/core_tree_metadata_complete.json').read_text());lookup={x['summary']['sample_id']:x for x in raw}
spec={'SRR10507650':('Aspergillus_floccosus','CBS 116.37','SAMN05878339'),'SRR8399223':('Aspergillus_aureoterreus','CBS 503.65','SAMN05518762'),'SRR35107300':('Aspergillus_hortae','IBT 26384 = CBS 124230','SAMN18260314')};anchors=[];flags=[]
for s,(species,strain,bio) in spec.items():
    b=next(x for x in lookup[s]['ncbi_biosamples'] if x['accession']==bio);attrs={x['attribute_name']:x['value'] for x in b['attributes']};assert attrs.get('type-material') and 'Genome Institute' in b['owner']
    anchors.append({'sample':s,'species':species,'culture':strain,'biosample':bio,'biosample_title':b['title'],'owner':b['owner'],'type_material':attrs['type-material'],'evidence':'Dedicated JGI culture-genome BioSample, matching collection strain and published ex-type identity. Genome-source/run mapping remains as documented in the assembly manifest.','culture_alias_source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7426331/'})
for d in raw:
    s=d['summary']['sample_id']
    if s in spec or not d['summary'].get('type_material'):continue
    flags.append({'sample':s,'strain':d['summary'].get('strain'),'isolate':d['summary'].get('isolate'),'type_material_claim':d['summary']['type_material'],'status':'not used as a type anchor','reason':'Clinical isolate collection metadata conflicts with treating the named reference culture/type-material field as verified strain identity; resolve provenance first.','biosamples':d['summary'].get('biosample_accessions')})
(P/'analysis/genome_type_anchor_audit.json').write_text(json.dumps({'accepted_genome_type_anchors':anchors,'withheld_type_claims':flags,'policy':'Do not infer type identity from an organism name, a type-material field alone, or marker proximity. No external reference records occur in reduced trees; three independently documented genome type cultures are present.'},indent=2));print('ANCHORS',len(anchors),'WITHHELD',len(flags))
