from pathlib import Path
import subprocess,os,time,json,hashlib,signal,datetime,sys
OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT))
from terreus_tree_helpers import parse_newick,support,strong
EXE='/home/benjamin/miniforge3/envs/klebsiella_recombination/bin/iqtree2'
PROTECTED='kleb-global-100gb-joint7830recovery1-20261118t003000z.service'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(name,d):
 p=OUT/name;t=p.with_suffix(p.suffix+'.pending');t.write_text(json.dumps(d,indent=2));t.replace(p)
def available():return int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
def protected():
 s=subprocess.check_output(['systemctl','--user','show',PROTECTED,'-p','MemoryMax','-p','ActiveState','-p','MemoryCurrent'],text=True)
 d=dict(x.split('=',1) for x in s.splitlines());assert d['MemoryMax']=='100000000000';return d
def main():
 cg=Path('/sys/fs/cgroup')/Path('/proc/self/cgroup').read_text().strip().split('::',1)[1].lstrip('/')
 limits={k:(cg/k).read_text().strip() for k in ['memory.max','memory.swap.max','cpu.max']}
 assert limits['memory.max']=='2147483648' and limits['memory.swap.max']=='0' and limits['cpu.max']=='100000 100000'
 assert available()>14*1024**3
 design=json.loads((OUT/'data/design.json').read_text());inp=OUT/'inputs/complete230.fasta';parts=OUT/'inputs/partitions.nex'
 for p,key in [(inp,'alignment_sha256'),(parts,'partitions_sha256')]:assert hashlib.sha256(p.read_bytes()).hexdigest()==design[key]
 version=subprocess.check_output([EXE,'--version'],text=True);assert 'version 2.4.0 ' in version
 prefix=OUT/'runs/complete230';argv=[EXE,'-s',str(inp),'-p',str(parts),'-m','MFP+MERGE','-B','1000','-alrt','1000','-bnni','-wbt','-seed','20261118','-nt','1','-o',design['outgroup'],'--prefix',str(prefix)]
 rec={'state':'running','started':now(),'argv':argv,'limits':limits,'protected_before':protected(),'version':version,'alignment_sha256':design['alignment_sha256']};save('status.json',rec)
 env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
 with (OUT/'runs/complete230.stdout').open('w') as f:
  child=subprocess.Popen(argv,stdout=f,stderr=subprocess.STDOUT,env=env,start_new_session=True)
  while child.poll() is None:
   if available()<12*1024**3:
    os.killpg(child.pid,signal.SIGTERM)
    try:child.wait(timeout=5)
    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
    raise RuntimeError('Stopped owned analysis because memory headroom fell below 12 GiB')
   time.sleep(2)
 assert child.returncode==0,('IQ-TREE exit',child.returncode)
 expected={l[1:].strip() for l in inp.read_text().splitlines() if l.startswith('>')}
 tree=parse_newick(Path(str(prefix)+'.treefile').read_text());assert {t.name for t in tree.leaves()}==expected
 n=0
 for line in Path(str(prefix)+'.ufboot').open():
  if line.strip():assert {t.name for t in parse_newick(line).leaves()}==expected;n+=1
 assert n==1000
 nodes=[]
 def visit(x):
  if x.children:nodes.append(x)
  for c in x.children:visit(c)
 visit(tree);labelled=[x for x in nodes if '/' in x.name];assert labelled
 report=Path(str(prefix)+'.iqtree').read_text();assert 'Date and time:' in report
 rec.update({'state':'complete','finished':now(),'tips':len(expected),'bootstrap_trees':n,'supported_internal_labels':len(labelled),'strong_internal_nodes':sum(strong(x) for x in nodes),'protected_after':protected(),'tree_sha256':hashlib.sha256(Path(str(prefix)+'.treefile').read_bytes()).hexdigest(),'memory_peak_bytes':int((cg/'memory.peak').read_text())})
 save('data/new_tree_validation.json',rec);save('status.json',rec);print('COMPLETE',len(expected),'tips',n,'bootstrap trees',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save('status.json',{'state':'failed','at':now(),'reason':str(e)});raise
