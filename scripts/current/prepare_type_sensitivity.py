from pathlib import Path
import json
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis';T=P/'tree_inputs'
occ={d['id']:d for d in json.loads((T/'locus_occupancy.json').read_text())};seq={};name=None
for line in (T/'combined_type_focused.fasta').read_text().splitlines():
    if line.startswith('>'):name=line[1:];seq[name]=''
    else:seq[name]+=line.strip()
exclude={k for k in seq if occ[k]['available_loci']<2 or k=='REF_453_Aspergillus_neoniveus_CBS26173'}
(T/'combined_type_focused_sensitivity.fasta').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in seq.items() if k not in exclude))
(T/'type_sensitivity_design.json').write_text(json.dumps({'input_tips':len(seq),'retained_tips':len(seq)-len(exclude),'excluded':sorted(exclude),'reason':'Reduce sparse-locus tips and remove composition-heterogeneous Neonivei outgroup; retain all nominal Terrei type-strain records and three other Flavipedes outgroups. Same trimmed sites, fixed partition models from full analysis.'},indent=2))
print('TYPE SENSITIVITY',len(seq)-len(exclude),'EXCLUDED',sorted(exclude))
