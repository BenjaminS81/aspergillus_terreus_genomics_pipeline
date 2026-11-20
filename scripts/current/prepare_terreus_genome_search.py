from pathlib import Path
import csv,json,collections,re
BASE=Path(__file__).resolve().parent;OUT=BASE/'taxonomy_reference_build/2026-11-10';P=OUT/'analysis'
local=json.loads((P/'local_metadata.json').read_text()); samples={d['sample'] for d in local.values()}
src=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\terreus_qc_check\complete_data\complete_genome_qc_summary.tsv')
rows=list(csv.DictReader(src.open(),delimiter='\t'))
mapping=[]
for sample in sorted(samples):
    matches=[d for d in rows if d['sample']==sample]
    if not matches:
        root=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
        paths=list((root/'results/0_2_repeat_mask_funannotate').glob('*/'+sample+'/contigs_renamed.fasta'))
        if not paths:paths=list((root/'results/0_1_spades').glob('*/'+sample+'/scaffolds.fasta'))
        if len(paths)!=1:raise ValueError((sample,len(matches),[str(p) for p in paths]))
        p=paths[0];d={'source_fasta':'/mnt/d/'+p.as_posix()[3:],'species_group':p.parent.parent.name,'data_class':'newer_source_not_in_prior_qc_table'}
    elif len(matches)==1:d=matches[0]
    else:raise ValueError((sample,len(matches)))
    winpath=Path('D:/'+d['source_fasta'].removeprefix('/mnt/d/'))
    if not winpath.is_file():
        root=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
        paths=list((root/'results/0_2_repeat_mask_funannotate').glob('*/'+sample+'/contigs_renamed.fasta'))
        if len(paths)!=1:raise ValueError(('missing_source',sample,str(winpath),[str(p) for p in paths]))
        d['source_fasta']='/mnt/d/'+paths[0].as_posix()[3:];d['data_class']+='; extant annotation input substituted for unavailable old QC path'
    mapping.append({'sample':sample,'subject':d['source_fasta'],'previous_species_group':d['species_group'],'qc_source':str(src),'data_class':d['data_class']})
(P/'genome_sources.json').write_text(json.dumps(mapping,indent=2))
refs=json.loads((P/'reference_metadata.json').read_text());r=next(d for d in refs if d['reference_id']=='REF_614_Aspergillus_terreus_CBS60165')
(P/'genome_search_queries.fasta').write_text(''.join('>'+l+'_'+r[l+'GB']+'\n'+r[l+'_sequence']+'\n' for l in ['BenA','CaM','RPB2']))
print('GENOME SOURCES',len(mapping))
