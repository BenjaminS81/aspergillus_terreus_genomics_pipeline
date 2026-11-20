from pathlib import Path
import json,subprocess,collections
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
calls=json.loads((P/'variant_genomic_markers/locus_calls.json').read_text());sources=json.loads((P/'variant_search_manifest.json').read_text());flag={d['sample'] for d in calls if d['status']=='no_strong_candidate_hold'}
O=P/'variant_rescue_blast';O.mkdir(exist_ok=True)
for d in sources:
    if d['sample'] not in flag:continue
    target=O/(d['sample']+'.tsv');cmd=['/usr/bin/blastn','-query',str(P/'rescue_type_queries.fasta'),'-subject',d['subject'],'-task','blastn','-evalue','1e-15','-dust','no','-perc_identity','70','-max_target_seqs','10000','-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen','-out',str(target)]
    r=subprocess.run(cmd,text=True,capture_output=True);(O/(d['sample']+'.log')).write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    d['output']=str(target.relative_to(P));d['rescue_command']=cmd;print('RESCUED VARIANT QUERY',d['sample'],d['sample_aliases'],flush=True)
(P/'variant_search_rescued_manifest.json').write_text(json.dumps(sources,indent=2))
