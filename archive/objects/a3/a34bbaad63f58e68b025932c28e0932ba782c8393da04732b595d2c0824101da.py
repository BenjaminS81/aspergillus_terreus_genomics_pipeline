from pathlib import Path
import sys,json,csv,hashlib,ast,math,struct,zlib,os
OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'comparative_visualization/2026-11-13/vendor'))
from circos_review_common import constrain_process
LIMITS=constrain_process()
from terreus_tree_helpers import parse_newick,support,prune
import pymupdf
import numpy as np
from PIL import Image
DPI=600;S=DPI/72;ROWS=128
src=ROOT/'stakeholder_release/2026-11-17_v3_taxonomy_review/scripts/tiff_writer_frozen.py'
nodes=[n for n in ast.parse(src.read_text()).body if isinstance(n,ast.FunctionDef)]
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(src),'exec'))
(OUT/'print_600dpi').mkdir(exist_ok=True)
records=[]
for p in sorted((OUT/'figures').glob('complete230*.pdf')):
 with pymupdf.open(p) as doc:
  assert len(doc)==1;page=doc[0];assert not page.get_image_info();assert all(doc.extract_font(f[0])[3] for f in page.get_fonts())
  r={'pdf':p.name,'sha256':sha(p),'pages':1,'all_vector':True,'fonts_embedded':True}
  r.update(write_tiff(page,OUT/'print_600dpi'/(p.stem+'.tif')));records.append(r)
  # Close-up checks at readable scale; full-page preview is produced by renderer.
  rect=page.rect;clip=pymupdf.Rect(0,0,rect.width,min(600,rect.height))
  page.get_pixmap(matrix=pymupdf.Matrix(1,1),clip=clip).save(OUT/'qa'/(p.stem+'_top.png'))
 print('600 dpi export',p.name,flush=True)
new=parse_newick((OUT/'runs/complete230.treefile').read_text());keep={t.name for t in new.leaves()};old=prune(parse_newick((OUT/'runs/full283_existing/combined_all_strains.treefile').read_text()),keep)
def splits(tree):
 result=set()
 def visit(n):
  a={n.name} if not n.children else set().union(*(visit(c) for c in n.children));b=keep-a
  if len(a)>1 and len(b)>1:result.add(min((tuple(sorted(a)),tuple(sorted(b))),key=lambda x:(len(x),x)))
  return a
 visit(tree);return result
a,b=splits(old),splits(new);comparison={'tips':len(keep),'RF':len(a^b),'denominator':len(a)+len(b),'normalised_RF':len(a^b)/(len(a)+len(b)),'scope':'new complete230 ML tree versus original full283 ML tree pruned to the identical 230 entries; sensitivity to cohort filtering and stochastic inference'}
(OUT/'data/full_vs_complete_RF.json').write_text(json.dumps(comparison,indent=2))
rows=[]
def visit(n):
 if n.children:
  x,y=support(n);rows.append({'node':len(rows)+1,'SH_aLRT':x,'UFBoot':y,'descendant_count':len(n.leaves()),'descendant_tips':';'.join(t.name for t in n.leaves())})
  for c in n.children:visit(c)
visit(new)
with (OUT/'data/branch_support.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
manifest={'new_figures':records,'full_tree':'Verified original IQ-TREE 2.4.0 run; original polished PDFs copied without modification','comparison':comparison}
log=(OUT/'runs/complete230.log').read_text()+'\n'+(OUT/'runs/complete230.iqtree').read_text();manifest['iqtree_diagnostics']=[l.strip() for l in log.splitlines() if any(t in l.lower() for t in ['warning','correlation','convergen','composition chi2','best partition model'])]
(OUT/'data/iqtree_diagnostics.json').write_text(json.dumps(manifest['iqtree_diagnostics'],indent=2))
calls=json.loads((ROOT/'taxonomy_reference_build/2026-11-10/analysis/genomic_markers/locus_calls.json').read_text())
selected=[{k:v for k,v in r.items() if k!='sequence'} for r in calls if 'LOCAL_'+r['sample'] in keep]
assert len(selected)==225 and all(r.get('strong_candidate_count',0)<=1 for r in selected)
(OUT/'data/selected_genome_marker_evidence.json').write_text(json.dumps(selected,indent=2))
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
(OUT/'README.txt').write_text('''THREE-GENE TREES — 2026-11-18 v1

A. Full reference-inclusive tree (283 tips)
Already inferred with IQ-TREE 2.4.0, MFP+MERGE, 1000 UFBoot, 1000 SH-aLRT,
-bnni, seed 20261110. All 1000 bootstrap trees and the exact plotted topology
were verified. Existing polished single-page PDFs and numerical outputs are
included under figures/full283_* and runs/full283_existing/. This was not rerun.

B. NEW complete-locus tree (230 tips)
75 genome-derived entries + 155 literature reference entries, including 29
reference type entries. Every entry represents BenA, CaM and RPB2. All eligible
entries retained; no downsampling and no substitution from other strains.
Same 2076 aligned sites: BenA 1-515, CaM 516-1062, RPB2 1063-2076.
Presence criterion >=120 ACGT per locus; actual minima are 392/434/637 bases.
Complete locus occupancy does not mean gap-free/full-length genes, proven
single-copy status, independent cultures or newly resolved orthology.
Existing genome search calls: 223 single strong candidates, two partial single
candidates (BenA in IFO6365 and SRR9691929); no selected call has multiple strong
candidates. Annotation mapping is incomplete, so this is not universal proof of
one-copy genomic orthology. See data/selected_genome_marker_evidence.json.

IQ-TREE 2.4.0: -p partitions.nex -m MFP+MERGE -B 1000 -alrt 1000 -bnni
-wbt -seed 20261118 -nt 1. Root display uses A. flavipes NRRL302, as before.
Exact arguments, versions, hashes and support validation: data/new_tree_validation.json.
RAM cap 2 GiB, one CPU, no swap, low priority. The 100 GB protected job is unchanged.

FIGURES
complete230_vertical.pdf: branch lengths in substitutions per site.
complete230_circular.pdf: cladogram; branch lengths not represented.
complete230_colour_legend.pdf: shared colours, type markers, support convention.
All are standalone single-page vector PDFs; SVG and native-size 600 dpi TIFF
exports supplied. QA PNGs are screen previews only.
Dots mark SH-aLRT >=80 and UFBoot >=95; all numerical values in Newick and
data/branch_support.tsv. Type-clade membership recalculated on the NEW tree;
T codes are not transferable to the earlier tree's C codes. Reviewed species
labels are inherited metadata, not newly asserted taxonomic determinations.

Bootstrap support is variable: 66 internal nodes pass both display thresholds;
IQ-TREE flags 78 near-zero internal branches and one composition-test failure.
Identical marker sequences are retained as separate labels; their within-group
ordering is not evidence of resolved strain relationships. Diagnostics are saved.

COVERAGE LIMITATIONS
Current A. jilinensis and A. ralucae reference type records lack all three loci
and are excluded. A. terreus var. subfloccosus has no complete reference entry.
HNGD-TM15 has BenA only and is excluded. SRR33623979 has no accepted marker calls.
No call is treated as proof of biological gene absence.
All selection and exclusion details: data/selection_audit.tsv and data/design.json.
This new tree does not expand the existing 66-entry core-genome backbone.

METHOD DOCUMENTATION
https://iqtree.github.io/doc/Command-Reference
''')
print('Package complete',comparison,flush=True)
