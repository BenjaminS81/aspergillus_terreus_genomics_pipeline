"""Validate GFF CDS coordinates/phases and translations before primary-protein export.
GFF3 phases follow transcript order; only the first CDS phase trims the joined CDS.
"""
from pathlib import Path
import sys,json,hashlib,collections,datetime,urllib.parse,warnings
from Bio.Seq import Seq
from Bio.Data import CodonTable
from worker import OUT,PREP,lp,load,save,sha
warnings.filterwarnings('error',category=Warning)
def fasta(p):
 name=None;parts=[]
 with Path(p).open() as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,''.join(parts).upper()
    name=line[1:].split()[0];parts=[]
   else:parts.append(line.strip())
 if name is not None:yield name,''.join(parts).upper()
def seqdict(p):
 d={}
 for h,s in fasta(p):
  if h in d:raise ValueError('Duplicate FASTA ID '+h)
  d[h]=s
 return d
def join_cds(parts,genome):
 if not parts:return '',[],['missing_CDS']
 errors=[];strands={x[3] for x in parts};contigs={x[0] for x in parts}
 if len(strands)!=1 or len(contigs)!=1 or not strands<={'+','-'}:return '',[],['inconsistent_CDS_location']
 strand=next(iter(strands));parts=sorted(parts,key=lambda x:x[1],reverse=strand=='-');chunks=[];length=0;prev=None
 for i,(contig,start,end,st,phase) in enumerate(parts):
  if contig not in genome or start<1 or end<start or end>len(genome.get(contig,'')):errors.append('invalid_CDS_coordinates');continue
  if phase not in {'0','1','2'}:errors.append('invalid_phase');continue
  phase=int(phase)
  if i and phase!=(-length)%3:errors.append('phase_discontinuity')
  if prev and (end>=prev[0] if strand=='-' else start<=prev[1]):errors.append('overlapping_CDS_segments')
  seq=genome[contig][start-1:end]
  if strand=='-':seq=str(Seq(seq).reverse_complement())
  if i==0:seq=seq[phase:]
  chunks.append(seq);length+=len(seq);prev=(start,end)
 return ''.join(chunks),parts,sorted(set(errors))
def selftest():
 # Split codon across two exons: phase 2 on second exon; don't delete its bases.
 s,p,e=join_cds([('c',1,4,'+','0'),('c',8,12,'+','2')],{'c':'ATGANNNAATAA'})
 assert s=='ATGAAATAA' and not e
 s,p,e=join_cds([('c',1,5,'-','2'),('c',9,12,'-','0')],{'c':'TTATTNNNTCAT'})
 assert s=='ATGAAATAA' and not e
 assert 'phase_discontinuity' in join_cds([('c',1,4,'+','0'),('c',8,12,'+','1')],{'c':'ATGANNNAATAA'})[2]
 save(OUT/'reports/CDS_validator_tests.json',{'positive_split_codon':'PASS','negative_split_codon':'PASS','wrong_phase_detected':'PASS'})
def validate(sample,a,proposal):
 genome_path=lp(a['scaffolds_file']);protein_path=lp(a['protein_file']);gff_path=lp(a['gff'])
 assert sha(protein_path)==a['file_sha256'] and sha(gff_path)==a['gff_sha256'],'Selected annotation changed'
 input_lock=load(PREP/'data/current_input_hashes.json')
 expected_genome=input_lock[str(a['scaffolds_file']).lower()]['file_sha256']
 assert sha(genome_path)==expected_genome,'Annotation genome changed'
 genome=seqdict(genome_path);proteins=seqdict(protein_path)
 cds=collections.defaultdict(list);parents={};gff_errors=[]
 with gff_path.open() as f:
  for lineno,line in enumerate(f,1):
   if line.startswith('#'):continue
   fields=line.rstrip().split('\t')
   if len(fields)!=9:continue
   contig,src,kind,start,end,score,strand,phase,attrs=fields
   if kind not in {'CDS','mRNA','transcript'}:continue
   attrs={k:urllib.parse.unquote(v) for k,v in (t.split('=',1) for t in attrs.split(';') if '=' in t)}
   if kind in {'mRNA','transcript'}:
    if 'ID' not in attrs or 'Parent' not in attrs or ',' in attrs['Parent']:gff_errors.append({'line':lineno,'issue':'missing_or_multiple_gene_parent'});continue
    if attrs['ID'] in parents:gff_errors.append({'line':lineno,'issue':'duplicate_transcript_ID'})
    parents[attrs['ID']]=attrs['Parent']
   else:
    for parent in attrs.get('Parent','').split(','):
     if parent:cds[parent].append((contig,int(start),int(end),strand,phase))
 orphan_cds=sorted(set(cds)-set(proteins))
 if orphan_cds:gff_errors.append({'issue':'CDS_without_protein','count':len(orphan_cds),'examples':orphan_cds[:20]})
 rows=[];valid=collections.defaultdict(list);sequences={};counts=collections.Counter()
 for pid,aa in proteins.items():
  errors=[];flags=[]
  dna,parts,errors=join_cds(cds.get(pid,[]),genome)
  if pid not in parents:errors.append('missing_gene_parent')
  if not dna:errors.append('empty_CDS')
  translated=str(Seq(dna[:len(dna)//3*3]).translate(table=1)) if dna else ''
  aa=aa.removesuffix('*');translated=translated.removesuffix('*')
  if len(dna)%3:flags.append('incomplete_terminal_codon')
  if dna[:3]!='ATG':flags.append('non_ATG_or_partial_start')
  if dna[-3:] not in {'TAA','TAG','TGA'}:flags.append('missing_terminal_stop')
  if parts and parts[0][4]!='0':flags.append('partial_initial_codon')
  match='exact' if translated==aa else 'mismatch'
  if match=='mismatch' and translated and aa and aa[0]=='M' and translated[1:]==aa[1:] and dna[:3] in CodonTable.unambiguous_dna_by_id[1].start_codons:
   match='valid_alternative_initiator';flags.append('alternative_initiator_M')
  if match=='mismatch':errors.append('CDS_protein_translation_mismatch')
  if '*' in aa or '*' in translated:errors.append('internal_stop')
  if set(aa)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):errors.append('invalid_amino_acid')
  if not aa:errors.append('empty_protein')
  errors=sorted(set(errors));row={'protein_id':pid,'gene_id':parents.get(pid),'aa_length':len(aa),'CDS_length':len(dna),'CDS_parts':len(parts),'translation_match':match,'errors':errors,'flags':flags}
  rows.append(row);counts['proteins']+=1;counts['valid']+=not errors;counts['partial_or_noncanonical']+=bool(flags)
  for err in errors:counts[err]+=1
  if not errors:valid[parents[pid]].append(row);sequences[pid]=aa
 # No primary export is allowed to masquerade as an intact release if models failed.
 selected=[]
 for gene,rr in sorted(valid.items()):
  rr.sort(key=lambda r:(-r['aa_length'],r['protein_id']));selected.append(rr[0]['protein_id'])
 folder=OUT/'runs/annotation_release'/sample;folder.mkdir(parents=True,exist_ok=True)
 with (folder/'primary_validated.faa').open('w') as f:
  for pid in selected:
   f.write('>'+pid+' gene='+parents[pid]+'\n');s=sequences[pid]
   for i in range(0,len(s),80):f.write(s[i:i+80]+'\n')
 selected_set=set(selected)
 for row in rows:row['primary_selected']=row['protein_id'] in selected_set
 result={'sample':sample,'source_protein':str(protein_path),'source_gff':str(gff_path),'source_genome':str(genome_path),'protein_sha256':sha(protein_path),'gff_sha256':sha(gff_path),'annotation_genome_sha256':sha(genome_path),'primary_sha256':sha(folder/'primary_validated.faa'),'counts':dict(counts),'primary_genes':len(selected),'gff_errors':gff_errors,'assembly_relationship':proposal['assembly_relationship'],'technical_status':'PASS' if counts['valid']==counts['proteins'] and not gff_errors else 'MODEL_EXCEPTIONS','biological_completeness':'not inferred; partial models and excluded contigs remain flagged','model_rows':rows}
 save(folder/'validation.json',result)
 return {k:v for k,v in result.items() if k!='model_rows'}
def main():
 selftest();proposal_file=OUT/'data/resolved_annotation_proposals.json';proposals=load(proposal_file if proposal_file.exists() else PREP/'data/annotation_release_proposals.json');ann={a['protein_file']:a for a in load(PREP/'data/annotation_candidates.json')}
 result=[]
 for p in proposals:
  sample=p['sample']
  if not p['preferred_candidate']:
   result.append({'sample':sample,'technical_status':'ASSEMBLY_ANNOTATION_LINK_UNRESOLVED'});continue
  report=OUT/'runs/annotation_release'/sample/'validation.json'
  if report.exists():
   r=load(report);assert sha(lp(p['preferred_candidate']))==r['protein_sha256'];r={k:v for k,v in r.items() if k!='model_rows'}
  else:r=validate(sample,ann[p['preferred_candidate']],p)
  result.append(r);save(OUT/'reports/annotation_release_summary.json',result)
  print(sample,r['technical_status'],r.get('counts'),flush=True)
 save(OUT/'reports/annotation_release_summary.json',result)
 save(OUT/'reports/annotation_acceptance.json',{'sample_labels':len(result),'status_counts':dict(collections.Counter(r['technical_status'] for r in result)),'gate_status':'technical model validation only; unresolved assemblies and biological alias/metadata gates remain explicit'})
if __name__=='__main__':main()
