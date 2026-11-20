"""Prepare one primary isoform per explicit GFF parent; never merge paralogous genes.
Run only after assembly/CDS release validation; selection does not close that gate.
"""
from pathlib import Path
import argparse, hashlib, json, collections
p=argparse.ArgumentParser()
p.add_argument('--protein',type=Path,required=True)
p.add_argument('--gff',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
if a.output.exists():raise SystemExit('Output exists; use a new versioned directory')
parents={}
for line in a.gff.read_text().splitlines():
 if line.startswith('#'):continue
 f=line.split('\t')
 if len(f)!=9 or f[2] not in {'mRNA','transcript'}:continue
 attrs=dict(t.split('=',1) for t in f[8].split(';') if '=' in t)
 if 'ID' in attrs and 'Parent' in attrs:
  if ',' in attrs['Parent'] or attrs['ID'] in parents:raise SystemExit('Ambiguous transcript parent')
  parents[attrs['ID']]=attrs['Parent']
def fasta(path):
 h=None;ss=[]
 for line in path.read_text().splitlines():
  if line.startswith('>'):
   if h is not None:yield h,''.join(ss).upper()
   h=line[1:];ss=[]
  else:ss.append(line.strip())
 if h is not None:yield h,''.join(ss).upper()
bygene=collections.defaultdict(list);seen=set();reject=[];total=0;dna=0
for h,s in fasta(a.protein):
 total+=1;pid=h.split()[0]
 if pid in seen:raise SystemExit('Duplicate protein ID: '+pid)
 seen.add(pid)
 if pid not in parents:raise SystemExit('Unmapped protein/GFF transcript ID: '+pid)
 gene=parents[pid];dna+=set(s)<=set('ACGTN')
 terminal_stop=s.endswith('*');s=s[:-1] if terminal_stop else s
 reason='empty' if not s else 'internal_stop' if '*' in s else 'invalid_amino_acid' if set(s)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO') else ''
 record={'protein_id':pid,'gene_id':gene,'length':len(s),'terminal_stop_removed':terminal_stop,'sequence':s}
 if reason:reject.append({**record,'reason':reason})
 else:bygene[gene].append(record)
if not total or dna/total>.9:raise SystemExit('Empty or predominantly nucleotide-only input; refuses amino-acid export')
selected=[];mapping=[]
for gene,rr in sorted(bygene.items()):
 rr.sort(key=lambda r:(-r['length'],r['protein_id']));selected.append(rr[0])
 for r in rr:mapping.append({k:v for k,v in {**r,'selected':r is rr[0]}.items() if k!='sequence'})
a.output.mkdir(parents=True)
with (a.output/'primary.faa').open('w') as f:
 for r in selected:
  f.write('>'+r['protein_id']+' gene='+r['gene_id']+'\n')
  for i in range(0,len(r['sequence']),80):f.write(r['sequence'][i:i+80]+'\n')
def sha(path):
 with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
(a.output/'selection.json').write_text(json.dumps({'input_protein_sha256':sha(a.protein),'input_gff_sha256':sha(a.gff),'primary_sha256':sha(a.output/'primary.faa'),'input_proteins':total,'selected_genes':len(selected),'mapping':mapping,'rejected':[{k:v for k,v in r.items() if k!='sequence'} for r in reject],'assembly_and_CDS_validation':'must be supplied separately; not established by this isoform export'},indent=2))
print(json.dumps({'input':total,'primary':len(selected),'rejected':len(reject),'output':str(a.output)}))
