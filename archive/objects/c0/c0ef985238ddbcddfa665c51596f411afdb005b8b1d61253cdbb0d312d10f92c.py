from pathlib import Path
import json,subprocess,collections,re
BASE=Path(__file__).resolve().parent;OUT=BASE/'taxonomy_reference_build/2026-11-10';P=OUT/'analysis'
ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
local=json.loads((P/'local_metadata.json').read_text());core={d['sample'] for d in local.values()}
paths=[]
# Enumerate assembly output locations without traversing millions of unrelated intermediates.
for pat in ['data/assemblies/*/*/*_genomic.fna','data/reference/*_genomic.fna','results*/0_1_spades/*/*/scaffolds.fasta','results*/0_1_spades/*/*/contigs.fasta','results*/0_2_repeat_mask_funannotate/*/*/contigs_renamed.fasta','results*/0_2_repeat_mask_funannotate/*/*/cleaned_masked_contigs.fasta','results*/0_2_repeat_mask_funannotate/*/*/*/annotate_results/*.fasta']:
    paths.extend(ROOT.glob(pat))
assets=[]
for p in paths:
    rel=p.relative_to(ROOT).as_posix();parts=rel.split('/');name=p.name;sample=None;priority=99
    if name=='cds_from_genomic.fna':continue
    if rel.startswith('data/assemblies/'):
        sample=parts[3];priority=0;kind='public_reference_genome'
    elif rel.startswith('data/reference/'):
        sample='NIH2624' if 'GCA_000149615.1' in name else p.stem;priority=2;kind='reference_copy'
    elif '/0_1_spades/' in rel:
        sample=parts[3];priority=3 if name=='scaffolds.fasta' else 4;kind='spades_assembly'
    elif '/0_2_repeat_mask_funannotate/' in rel:
        sample=parts[3];priority=1 if name=='contigs_renamed.fasta' else 5 if len(parts)==5 else 7;kind='annotation_genome_input'
    else:kind='unmapped_assembly_candidate'
    if sample=='ATCC_20541':sample='ATCC20541'
    if sample and sample.casefold() in {c.casefold() for c in core}:sample=next(c for c in core if c.casefold()==sample.casefold())
    assets.append({'sample':sample,'windows_path':str(p),'subject':'/mnt/d/'+p.as_posix()[3:],'relative_path':rel,'bytes':p.stat().st_size,'priority':priority,'source_kind':kind})
mapping=[]
for sample in sorted({d['sample'] for d in assets if d['sample']}):
    options=[d for d in assets if d['sample']==sample and d['bytes']>0]
    if not options:continue
    d=sorted(options,key=lambda x:(x['priority'],x['relative_path']))[0].copy();d['in_core_tree']=sample in core
    mapping.append(d)
missing=core-{d['sample'] for d in mapping};assert not missing,missing
for d in assets:d['selected_canonical']=any(d['subject']==m['subject'] for m in mapping)
(P/'assembly_asset_inventory.json').write_text(json.dumps(assets,indent=2))
(P/'genome_sources.json').write_text(json.dumps(mapping,indent=2))
refs=json.loads((P/'reference_metadata.json').read_text());r=next(d for d in refs if d['reference_id']=='REF_614_Aspergillus_terreus_CBS60165')
(P/'genome_search_queries.fasta').write_text(''.join('>'+l+'_'+r[l+'GB']+'\n'+r[l+'_sequence']+'\n' for l in ['BenA','CaM','RPB2']))
print('ASSETS',len(assets),'EMPTY',sum(d['bytes']==0 for d in assets),'GENOMES',len(mapping),'OUTSIDE_CORE',[d['sample'] for d in mapping if not d['in_core_tree']],'UNMAPPED',[d['relative_path'] for d in assets if not d['sample']])
