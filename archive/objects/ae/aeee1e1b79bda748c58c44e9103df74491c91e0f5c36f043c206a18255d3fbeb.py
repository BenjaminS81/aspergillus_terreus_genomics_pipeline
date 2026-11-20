from pathlib import Path
import sys,subprocess,json
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from circos_review_common import constrain_process
limits=constrain_process()
node=r'C:\Users\q287js\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe'
script=Path(__file__).with_name('build_workbook.mjs')
rc=subprocess.call([node,'--max-old-space-size=768',str(script)])
(script.parent.parent/'qa/workbook_resource_limits.json').write_text(json.dumps(limits,indent=2))
sys.exit(rc)
