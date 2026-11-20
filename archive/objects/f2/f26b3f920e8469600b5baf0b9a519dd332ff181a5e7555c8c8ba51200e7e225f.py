"""Pin current core alignment and prepare, without executing, the support command."""
from pathlib import Path
import sys,json,hashlib
OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
constrain_process()
p=Path('D:/Illumina/David_Eisele/terreus_genotyping/results/0_4_mafft_iqtree/Sep28/concatenated.fa')
e=json.loads((ROOT/'review_reports/2026-11-10/audit_evidence.json').read_text())
rows=[];name=None;size=0
with p.open() as f:
 for line in f:
  if line.startswith('>'):
   if name is not None:rows.append({'id':name,'length':size})
   name=line[1:].split()[0];size=0
  else:size+=len(line.strip())
if name is not None:rows.append({'id':name,'length':size})
assert len(rows)==66 and {r['id'] for r in rows}==set(e['trees']['genome']['tips'])
assert all(r['length']==9012257 for r in rows)
with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
linux='/mnt/d/Illumina/David_Eisele/terreus_genotyping/results/0_4_mafft_iqtree/Sep28/concatenated.fa'
prefix='/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics/gate_preparation/2026-11-18_v1/runs/core_support/version_sensitivity_2.4.0'
record={'status':'prepared, not executed; version sensitivity unless original 3.1.1 is recovered','input_windows':str(p),'input_linux':linux,'input_sha256':digest,'shape':{'tips':66,'columns':9012257},'historical_exact_block_validation':'review_reports/2026-11-10/concatenation_order_verification.json: zero mismatches; current check confirms shape/tips, not a repeated full block reconstruction','argv':['/home/benjamin/miniforge3/envs/klebsiella_recombination/bin/iqtree2','-s',linux,'-m','GTR+F+R10','-B','1000','-alrt','1000','-bnni','-wbt','-seed','20261118','-nt','1','-mem','4G','--prefix',prefix],'prerequisites':['dedicated resource window; bounded cgroup; output directory absent/new','recheck input hash before launch','explicitly retain 2.4.0 version sensitivity designation or pin original 3.1.1','root selection remains separate; no -o assumption'],'followups':['supported partitioned/model sensitivity with verified recovered coordinates','composition failure investigation','amino-acid core after validated G4 proteins','RF/NMDS comparisons on identical 66 tips; no taxonomic authority ranking'],'original_log_evidence':e['trees']['genome']['log_evidence']}
(OUT/'data/core_support_job_plan.json').write_text(json.dumps(record,indent=2))
print('Core alignment pinned:',digest,'66 x 9012257; job not launched')
