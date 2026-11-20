from pathlib import Path
import json,csv,re,collections
ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping');P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis/orthology_crosscheck'
ogs=json.loads((P/'relevant_orthogroup_rows.json').read_text());mapped=json.loads((P/'validated_marker_gene_mapping.json').read_text());calls=json.loads((P.parent/'genomic_markers/locus_calls.json').read_text())
samples=[k.removesuffix('_transcripts') for k in next(iter(ogs.values())) if k!='Orthogroup'];out=[]
for sample in samples:
    folders=list((ROOT/'results/0_2_repeat_mask_funannotate').glob('*/'+sample+'/*/annotate_results'));assert folders,sample
    gff=next(sorted(folders)[0].glob('*.gff3'));genes={};mrna={}
    for line in gff.read_text().splitlines():
        if line.startswith('#'):continue
        a=line.split('\t')
        if len(a)!=9 or a[2] not in ['gene','mRNA']:continue
        attr=dict(x.split('=',1) for x in a[8].split(';') if '=' in x)
        if a[2]=='gene':genes[attr['ID']]={'contig':a[0],'start':int(a[3]),'end':int(a[4]),'strand':a[6],'name':attr.get('Name','')}
        else:mrna[attr.get('Parent','')]=attr
    for og,row in ogs.items():
        for member in [x.strip() for x in row[sample+'_transcripts'].split(',') if x.strip()]:
            gid=re.sub('-T\d+$','',member);g=genes.get(gid,{});m=mrna.get(gid,{})
            loci=[d['locus'] for d in mapped if d['sample']==sample and any(x['gene_id']==gid for x in d.get('genes',[]))]
            out.append({'sample':sample,'orthogroup':og,'transcript':member,'gene_id':gid,**g,'product':m.get('product',''),'dbxref':m.get('Dbxref',''),'validated_marker_overlap':','.join(loci),'tubulin_domain':any(x in m.get('Dbxref','') for x in ['PF00091','PF03953']),'RRM_domain':'PF00076' in m.get('Dbxref',''),'metallophosphoesterase_domain':'PF00149' in m.get('Dbxref',''),'annotation_source':str(gff)})
cols=list(dict.fromkeys(k for d in out for k in d));
with (P/'all_relevant_orthogroup_members.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(out)
(P/'all_relevant_orthogroup_members.json').write_text(json.dumps(out,indent=2))
summary=[]
for sample in samples:
    d={'sample':sample}
    for og in ogs:d[og+'_gene_models']=len({x['gene_id'] for x in out if x['sample']==sample and x['orthogroup']==og})
    for l in ['BenA','CaM','RPB2']:
        call=next(c for c in calls if c['sample']==sample and c['locus']==l);d[l+'_genomic_candidates']=call['strong_candidate_count'];d[l+'_annotation_models_matching_marker']=len(next(x for x in mapped if x['sample']==sample and x['locus']==l).get('genes',[]))
    summary.append(d)
(P/'gene_models_vs_genomic_loci.json').write_text(json.dumps(summary,indent=2))
print('ORTHOGROUP MEMBERS',len(out));print('NIH2624',json.dumps([d for d in out if d['sample']=='NIH2624'],indent=2))
