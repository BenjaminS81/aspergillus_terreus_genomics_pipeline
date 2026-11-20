"""Recover exact legacy DNA inputs and IDs from OrthoFinder's retained workspace."""
from pathlib import Path
import json,collections,re,warnings
from Bio.Seq import Seq
from worker import OUT,load,save,sha
from validate_annotations import fasta,seqdict
SRC=Path('/mnt/d/Illumina/David_Eisele/terreus_genotyping/results/0_3_orthofinder/terrei3/Results_Sep28')
def main():
 w=SRC/'WorkingDirectory';log=(SRC/'Log.txt').read_text();assert ' -d ' in log and '3.0.1b1' in log
 species=dict(x.split(': ',1) for x in (w/'SpeciesIDs.txt').read_text().splitlines() if ': ' in x)
 ids=dict(x.split(': ',1) for x in (w/'SequenceIDs.txt').read_text().splitlines() if ': ' in x)
 folder=OUT/'runs/recovered_orthofinder_DNA';folder.mkdir(exist_ok=True);records=[]
 for key,name in species.items():
  sample=name.removesuffix('_transcripts.fa');source=w/('Species'+key+'.fa');target=folder/name
  primary_path=OUT/'runs/annotation_release'/sample/'primary_validated.faa';primary=seqdict(primary_path) if primary_path.exists() else {}
  counts=collections.Counter();seen=set()
  with target.open('w') as f:
   for sid,seq in fasta(source):
    original=ids[sid];pid=original.split()[0];assert pid not in seen;seen.add(pid)
    counts['sequences']+=1;counts['bases']+=len(seq);counts['non_DNA_records']+=bool(set(seq)-set('ACGTURYSWKMBDHVN'))
    f.write('>'+original+'\n')
    for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')
    if pid in primary:
     counts['IDs_in_current_primary']+=1
     tr=str(Seq(seq[:len(seq)//3*3]).translate()).removesuffix('*')
     counts['direct_coding_translation_exact_primary_matches']+=tr==primary[pid]
  records.append({'sample':sample,'species_index':key,'source':str(source),'source_sha256':sha(source),'recovered_input':str(target),'recovered_sha256':sha(target),'counts':dict(counts),'comparison_policy':'Direct full-transcript coding translation is a conservative identity check, not gene-model reconstruction; UTR/phase differences can cause nonmatches.'})
  print('Recovered legacy DNA input',sample,counts['sequences'],flush=True)
 assert len(records)==66 and all(r['counts']['non_DNA_records']==0 for r in records)
 save(OUT/'reports/orthofinder_input_recovery.json',{'tool':'OrthoFinder 3.0.1b1','input_mode':'DNA (-d)','search_evidence':'blast_nuclDBSpecies* retained; original command includes -d','original_command':next(x for x in log.splitlines() if x.startswith('Command Line:')),'species_ids_sha256':sha(w/'SpeciesIDs.txt'),'sequence_ids_sha256':sha(w/'SequenceIDs.txt'),'log_sha256':sha(SRC/'Log.txt'),'input_files':records,'orthogroup_status':'Exact DNA inputs recovered; do not relabel as a protein-input analysis. Final species-tree/HOG completion and primary-protein sensitivity remain unverified.'})
if __name__=='__main__':main()
