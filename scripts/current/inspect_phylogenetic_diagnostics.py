from pathlib import Path
import json,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis';out={}
occupancy={d['id']:d for d in json.loads((P/'tree_inputs/locus_occupancy.json').read_text())}
for name in ['combined_all_strains','type_focused_sensitivity']:
    lines=(P/'trees'/(name+'.mldist')).read_text().splitlines();rows=[ln.split() for ln in lines[1:] if ln.strip()];names=[r[0] for r in rows];dist=[[float(x) for x in r[1:]] for r in rows];high=[]
    for i in range(len(rows)):
        for j in range(i):
            if dist[i][j]>=1:high.append({'tip1':names[i],'tip2':names[j],'distance':dist[i][j]})
    complete=[d for d in high if occupancy[d['tip1']]['available_loci']==occupancy[d['tip2']]['available_loci']==3]
    out[name]={'largest_distance':max(max(r) for r in dist),'pairs_distance_ge_1':high,'high_distance_pairs_with_both_tips_three_loci':len(complete)};print(name,'MAX',out[name]['largest_distance'],'HIGH',len(high),'BOTH COMPLETE',len(complete),'EXAMPLES',complete[:3])
(P/'pairwise_distance_diagnostics.json').write_text(json.dumps(out,indent=2))
tax=json.loads((P/'sample_taxonomic_review.json').read_text());print('SENSITIVITY',collections.Counter(d['full_vs_sensitivity'] for d in tax))
