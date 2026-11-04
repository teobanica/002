"""Complete, checkpointed A/B fourth-row scans for the 32 uniform AAA components."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,subprocess,os,hashlib
R=Path(__file__).resolve().parents[1];W=R/'stage_uniform'
for sub in ['inputs','outputs','checks']:(W/sub).mkdir(parents=True,exist_ok=True)
cs=json.loads((R.parent/'n8_aaa_complete/inputs/components.json').read_text())
cs=[dict(c,component_index=i) for i,c in enumerate(cs) if c['class_id'] not in (11,19,21,22)]
assert len(cs)==32
(W/'inputs/components.json').write_text(json.dumps(cs,separators=(',',':'))+'\n')
inp=W/'inputs/components.txt'
with inp.open('w') as f:
 f.write(str(len(cs))+'\n')
 for c in cs:
  f.write(f"{c['component_index']} {c['branch']} {c['dim']}\n")
  for row in c['rows']:f.write(' '.join(str(x) for e in row for x in e)+'\n')
exe=R/'checks/full_rows';signature=hashlib.sha256(exe.read_bytes()+inp.read_bytes()).hexdigest()
def run(c):
 i=c['component_index'];out=W/f'outputs/full_rows_{i:02d}.json';sig=W/f'checks/full_rows_{i:02d}.sha256';log=W/f'checks/full_rows_{i:02d}.log';tmp=out.with_suffix('.pending')
 if out.exists() and sig.exists() and sig.read_text().strip()==signature:return
 with log.open('w') as err:r=subprocess.run([str(exe),str(inp),str(i)],stdout=subprocess.PIPE,stderr=err)
 assert r.returncode==0,(i,r.returncode)
 d=json.loads(r.stdout);assert len(d)==1 and d[0]['id']==i
 tmp.write_bytes(r.stdout)
 os.replace(tmp,out);sig.write_text(signature+'\n')
 print('component',i,'class',c['class_id'],'rows',len(d[0]['families']),flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,cs))
print('All 32 uniform AAA components scanned.',flush=True)
