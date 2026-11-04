from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from itertools import permutations
import json,subprocess,os,hashlib
R=Path(__file__).resolve().parents[1];W=R/'stage_bbb'
doc=json.loads((W/'inputs/components.json').read_text());MOD=doc['modulus'];components=[];seen={};maps=[]
def canonical(H):
 best=None;cert=None
 for rp in permutations(range(3)):
  for p in range(8):
   cols=[tuple((H[r][j]-H[r][p]-H[rp[0]][j]+H[rp[0]][p])%MOD for r in rp) for j in range(8)]
   cp=sorted(range(8),key=lambda j:cols[j]);key=tuple(cols[j] for j in cp)
   if best is None or key<best:best=key;cert=dict(rows=rp,pivot_column=p,columns=cp)
 return best,cert
for s in doc['components']:
 if any('A' in p for p in s['profiles']):continue
 if s['dim']==0:
  key,cert=canonical([[e[0] for e in row] for row in s['rows']])
  if key in seen:maps.append(dict(source=s['id'],target=seen[key],transform=cert));continue
  seen[key]=dict(id=len(components),transform=cert)
 components.append(dict(s,id=len(components),source_id=s['id'],component_index=len(components)))
(W/'inputs/reduced_components.json').write_text(json.dumps(dict(modulus=MOD,components=components,point_equivalences=maps),separators=(',',':'))+'\n')
inp=W/'inputs/components.txt'
with inp.open('w') as f:
 f.write(str(len(components))+'\n')
 for c in components:
  f.write(f"{c['id']} 0 {c['dim']}\n")
  for row in c['rows']:f.write(' '.join(str(x) for e in row for x in e)+'\n')
exe=R/'checks/full_rows_general';signature=hashlib.sha256(exe.read_bytes()+inp.read_bytes()).hexdigest()
print('Scanning',len(components),'BBB components, including',len(seen),'distinct fixed triples.',flush=True)
def run(c):
 i=c['id'];out=W/f'outputs/full_rows_{i:02d}.json';sig=W/f'checks/full_rows_{i:02d}.sha256';log=W/f'checks/full_rows_{i:02d}.log'
 if out.exists() and sig.exists() and sig.read_text().strip()==signature:return
 with log.open('w') as err:r=subprocess.run([str(exe),str(inp),str(i),str(MOD)],stdout=subprocess.PIPE,stderr=err)
 assert r.returncode==0,(i,r.returncode)
 d=json.loads(r.stdout);assert len(d)==1 and d[0]['id']==i
 tmp=out.with_suffix('.pending');tmp.write_bytes(r.stdout);os.replace(tmp,out);sig.write_text(signature+'\n')
 print(i,'dim',c['dim'],'fourth row records',len(d[0]['families']),flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,components))
print('All BBB fourth-row scans finished.',flush=True)
