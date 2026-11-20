from pathlib import Path
import json,re,collections
ROOT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\Houbraken_genes');P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis/orthology_crosscheck'
def fasta(p):
    n=None;s=[]
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            if n:yield n,''.join(s).upper().replace('-','')
            n=line[1:].split()[0];s=[]
        else:s.append(line.strip())
    if n:yield n,''.join(s).upper().replace('-','')
raw=list(fasta(ROOT/'BENA56.fa'));meta=json.loads((P/'old_marker_genes.json').read_text());out=[]
for name,seq in fasta(ROOT/'merged/individual/BENA56_aligned.fa'):
    if not set(seq)&set('ACGT'):continue
    exact=[h for h,s in raw if s==seq and h.split('|')[0]==name]
    sample=re.sub(r'^(terreus|hortae|pseudoterreus|alabamensis|aureoterreus|floccosus)_','',name)
    genes=[x for x in meta if x['sample']==sample and x['locus']=='BenA' and x['original_header'] in exact]
    out.append({'aligned_label':name,'sample':sample,'exact_raw_headers':exact,'orthogroups':sorted({x['orthogroup'] for x in genes if x.get('orthogroup')}),'sequence_length':len(seq)})
(P/'actual_old_BenA_alignment_mapping.json').write_text(json.dumps(out,indent=2))
print('OLD ALIGNMENT',len(out),'MATCHED',sum(bool(x['exact_raw_headers']) for x in out),'OG',collections.Counter(str(x['orthogroups']) for x in out))
print('RAW EXTRA',[(h,len(s)) for h,s in raw if set(s)&set('ACGT') and not any(h in x['exact_raw_headers'] for x in out)])
