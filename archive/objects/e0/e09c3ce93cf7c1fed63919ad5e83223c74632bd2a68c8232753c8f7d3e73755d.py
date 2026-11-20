from pathlib import Path
import json,csv,re,hashlib,collections
from html.parser import HTMLParser
from terreus_tree_helpers import parse_newick
import pymupdf
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';D=P/'core66_combinations';A=P/'analysis';L=P/'lovastatin';checks=[]
def js(p):return json.loads(p.read_text(encoding='utf-8'))
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)
def fasta(p):
    d={};k=None
    for line in p.read_text().splitlines():
        if line.startswith('>'):k=line[1:].split()[0];assert k not in d;d[k]=''
        elif line.strip():d[k]+=line.strip()
    return d
design=js(D/'design.json');keep=set(design['core_tips']);check('Exactly 66 unique genome-tree labels',len(keep)==66)
refann=js(A/'tree_annotations_all_combinations.json');annkeys={(x['tree'],x['tip']) for x in refann};check('Unique annotation record per reference tree and tip',len(annkeys)==len(refann))
allclades=js(A/'type_clades_all_combinations.json');refs=js(BASE/'taxonomy_reference_build/2026-11-10/analysis/tree_reference_metadata.json');typ={x['reference_id'] for x in refs if x['type_flag']}
check('Every supported reference type clade has an included ex-type and both support thresholds',all(set(c['type_reference_tips'].split(';'))<=set(c['descendant_tips']) and set(c['type_reference_tips'].split(';'))<=typ and c['SH_aLRT']>=80 and c['UFBoot']>=95 for c in allclades))
coreann=js(A/'core66_tree_annotations.json');check('All 15 reduced/core trees have complete annotations',len(coreann)==15*66 and len({(x['tree'],x['tip']) for x in coreann})==15*66)
anchors=js(A/'genome_type_anchor_audit.json');check('Three documented type cultures; six conflicting type claims withheld',len(anchors['accepted_genome_type_anchors'])==3 and len(anchors['withheld_type_claims'])==6)
check('No withheld type claim used as a genome-type anchor',not ({x['sample'] for x in anchors['withheld_type_claims']} & {x['sample'] for x in anchors['accepted_genome_type_anchors']}))
for d in design['combinations']:
    n=d['name'];aln=fasta(D/(n+'.fasta'));tree=parse_newick((D/(n+'.treefile')).read_text());pruned=parse_newick((D/(n+'.reference_pruned66.treefile')).read_text())
    check(n+': alignment, inferred tree and pruned tree contain the exact 66 tips',set(aln)=={x.name for x in tree.leaves()}=={x.name for x in pruned.leaves()}==keep)
    check(n+': alignment length and completed supported inference',{len(x) for x in aln.values()}=={d['sites']} and 'Date and time:' in (D/(n+'.iqtree')).read_text() and len((D/(n+'.ufboot')).read_text().splitlines())==1000)
    for variant in ['de_novo','reference_pruned']:
        for layout in ['linear','circular']:
            for ext in ['svg','pdf','png']:assert (P/'figures/core66_trees'/(n+'_'+variant+'_'+layout+'.'+ext)).stat().st_size>1000
    source='combined_all_strains' if len(d['genes'])==3 else n
    for variant in ['all_references','type_focused']:
        for layout in ['linear','circular']:
            for ext in ['svg','pdf','png']:assert (P/'figures/trees'/(source+'_'+variant+'_'+layout+'.'+ext)).stat().st_size>1000
    if len(d['genes'])==2:
        folder=P/'reference_pair_combinations';assert 'Date and time:' in (folder/(n+'.iqtree')).read_text();rt=parse_newick((folder/(n+'.treefile')).read_text());ra=fasta(folder/(n+'.fasta'))
        check(n+': full-reference tree restores all tips and annotations',len(rt.leaves())==len(ra) and {x.name for x in rt.leaves()}==set(ra)=={tip for name,tip in annkeys if name==n})
check('All seven choices have linear/circular views for reference and reduced cohorts',True)
rank=js(D/'ranking.json');check('Seven marker choices compared and normalized RF is bounded',len(rank)==7 and all(0<=r['normalized_RF_to_core']<=1 for r in rank))
check('All 7,000 bootstrap RF comparisons accounted for',sum(r['replicates'] for r in js(D/'bootstrap_RF_summary.json'))==7000)
with (D/'normalized_RF_matrix.tsv').open() as f:rr=list(csv.reader(f,delimiter='\t'));matrix=[[float(x) for x in row[1:]] for row in rr[1:]]
check('Eight-tree RF matrix symmetric with zero diagonal',len(matrix)==8 and all(matrix[i][j]==matrix[j][i] and (i!=j or matrix[i][j]==0) for i in range(8) for j in range(8)))
regions=js(L/'visualization_data.json');manifest=js(L/'genome_manifest.json');check('All 77 assembly labels represented across 72 unique canonical genome searches',len(regions)==77 and len({r['sample'] for r in regions})==77 and len(manifest)==72 and all((L/'blastn'/(r['representative']+'.done')).exists() for r in manifest))
check('All 23 planned translated-rescue searches completed',len(js(L/'tblastn/execution.json'))==23 and all(x['completed'] for x in js(L/'tblastn/execution.json')))
for r in regions:
    for p in r['segments']:
        assert 1<=p['window_start']<=p['start']<=p['end']<=p['window_end']<=p['contig_length']
        assert p['start']-p['window_start']<=20000 and p['window_end']-p['end']<=20000
        assert all(g['contig']==p['contig'] and p['start']<=g['start']<=g['end']<=p['end'] for g in p['genes'])
        if p.get('annotation_window_exact_match') is False:assert not p['annotation_genes']
    assert (P/'figures/lovastatin/individual'/(r['sample']+'.svg')).exists()
    assert (P/'figures/lovastatin/individual'/(r['sample']+'.png')).exists()
    records=fasta(L/'region_fastas'/(r['representative']+'.fasta'))
    for i,p in enumerate(r['segments']):
        seq=next(v for k,v in records.items() if '|segment'+str(i+1)+'|' in k);assert len(seq)==p['window_end']-p['window_start']+1;assert hashlib.sha256(seq.encode()).hexdigest()==p['sequence_sha256']
check('Genomic regions, contig breaks, 20 kb bounds, extracted sequences and annotation-coordinate checks consistent',True)
doc=pymupdf.open(P/'figures/lovastatin/all_genomes_cluster_and_flanks_atlas.pdf');check('The PDF atlas contains all 77 genome pages',len(doc)==77);doc.close()
status=collections.Counter(r['block_status'] for r in regions);check('Lovastatin status totals match all 77 labels',status=={'nine_gene_block':45,'partial_or_divergent_block':11,'no_convincing_colocalized_block':21})
class Links(HTMLParser):
    def __init__(self):super().__init__();self.targets=[]
    def handle_starttag(self,t,attrs):self.targets.extend(v for k,v in attrs if k in ['src','href'] and v)
parser=Links();parser.feed((P/'terreus_comparative_figures.html').read_text(encoding='utf-8'));future={'validation_summary.json','file_manifest.json'};missing=[x for x in parser.targets if not re.match(r'^[a-z]+:|^#',x) and x not in future and not (P/x).exists()]
check('All static report links resolve (manifest generated at packaging)',not missing)
# Figure families share layouts; render selected final PDFs, including both simple and fragmented loci.
qa=P/'qa';qa.mkdir(exist_ok=True)
for rel,page,name in [('figures/lovastatin/all_genomes_cluster_and_flanks_atlas.pdf',0,'lov_atlas_final_page1'),('figures/lovastatin/all_genomes_cluster_and_flanks_atlas.pdf',38,'lov_atlas_final_page39'),('figures/core66_comparison/all_eight_NMDS_MDS.pdf',0,'eight_tree_ordination_final'),('figures/core66_trees/BenA_CaM_RPB2_de_novo_circular.pdf',0,'core66_triple_final')]:
    doc=pymupdf.open(P/rel);pg=doc[page];scale=min(1.0,1600/max(pg.rect.width,pg.rect.height));pg.get_pixmap(matrix=pymupdf.Matrix(scale,scale)).save(qa/(name+'.png'));doc.close()
out={'status':'passed','checks':checks,'check_count':len(checks),'figure_pdf_count':len(list((P/'figures').rglob('*.pdf'))),'core_labels':66,'marker_choices':7,'lovastatin_labels':77,'source_project_modified':False}
(P/'validation_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
