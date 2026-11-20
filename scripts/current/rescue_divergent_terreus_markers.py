from pathlib import Path
import json,subprocess,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
calls=json.loads((P/'genomic_markers/locus_calls.json').read_text());refs=json.loads((P/'validated_reference_metadata.json').read_text());sources=json.loads((P/'genome_search_manifest.json').read_text())
samples=sorted({c['sample'] for c in calls if c['status']!='single_strong_candidate'})
fa=P/'rescue_type_queries.fasta';fa.write_text(''.join('>'+l+'__'+d['reference_id']+'\n'+d[l+'_analysis_sequence']+'\n' for d in refs if d['type_flag'] or d['accepted_type_name'] for l in ['BenA','CaM','RPB2'] if d[l+'_analysis_sequence']))
out=P/'rescue_blast';out.mkdir(exist_ok=True)
manifest=[]
for sample in samples:
    d=next(d for d in sources if d['sample']==sample);target=out/(sample+'.tsv')
    cmd=['/usr/bin/blastn','-query',str(fa),'-subject',d['subject'],'-task','blastn','-evalue','1e-15','-dust','no','-perc_identity','70','-max_target_seqs','10000','-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen','-out',str(target)]
    r=subprocess.run(cmd,text=True,capture_output=True);(out/(sample+'.log')).write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    e=d.copy();e['output']=str(target.relative_to(P));e['rescue_command']=cmd;manifest.append(e);print('BROAD REFERENCE RESCUE',sample,flush=True)
original={d['sample']:d for d in sources};original.update({d['sample']:d for d in manifest})
(P/'genome_search_rescued_manifest.json').write_text(json.dumps(list(original.values()),indent=2))
