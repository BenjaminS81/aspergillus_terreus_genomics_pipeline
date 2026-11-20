from pathlib import Path
import json
ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping');P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
records=[]
for line in (ROOT/'data/assemblies/annotated/NIH2624/genomic.gff').read_text().splitlines():
    if line.startswith('#'):continue
    a=line.split('\t')
    if len(a)!=9 or a[2]!='mRNA':continue
    if (a[0]=='CH476594.1' and int(a[3])<=822461<=int(a[4])) or (a[0]=='CH476601.1' and int(a[3])<=1682843<=int(a[4])):
        records.append(a);print('\t'.join(a))
(P/'NIH2624_tubulin_NCBI_annotation.json').write_text(json.dumps(records,indent=2))
