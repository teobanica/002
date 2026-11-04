from pathlib import Path
import json,sys
from collections import Counter
R=Path(__file__).resolve().parents[1]
name=sys.argv[1] if len(sys.argv)>1 else 'all'
ds=json.loads((R/f'outputs/full_rows_{name}.json').read_text())
def implies(a,v):
 c,t,*free=v
 if any(free):return False
 if a['eq']:
  k,b=a['eq'][0]
  if t%k:return False
  c-=t//k*b
 elif t:return False
 return c%30==0
def contains(a,b):
 # Returns a map proving the row/locus a is contained in b.
 if a['free']>b['free']:return None
 if any(not implies(a,[c,k]) for k,c in b['eq']):return None
 w=2+a['free'];values=[]
 for k in range(b['free']):
  j=next(j for j,e in enumerate(b['row']) if e[2+k])
  assert b['row'][j][2+k]==1 and sum(b['row'][j][2:])==1
  values.append([a['row'][j][0]-b['row'][j][0],a['row'][j][1]-b['row'][j][1]]+a['row'][j][2:])
 for x,y in zip(a['row'],b['row']):
  z=y[:2]+[0]*a['free']
  for k,v in enumerate(values):
   for t in range(w):z[t]+=y[2+k]*v[t]
  if not implies(a,[u-v for u,v in zip(z,x)]):return None
 return values
out=[]
for d in ds:
 fs=d['families'];keep=[];maps=[]
 order=sorted(range(len(fs)),key=lambda i:(-(1-len(fs[i]['eq'])+fs[i]['free']),-fs[i]['free'],-(fs[i]['eq'][0][0] if fs[i]['eq'] else 0),i))
 for i in order:
  for j in keep:
   v=contains(fs[i],fs[j])
   if v is not None:maps.append({'source':i,'target':j,'substitution':v});break
  else:keep.append(i)
 out.append(dict(d,keep=keep,containment_maps=maps))
 print(d['id'],len(fs),'->',len(keep),dict(Counter((len(fs[i]['eq']),fs[i]['free']) for i in keep)),flush=True)
(R/f'outputs/reduced_rows_{name}.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
with (R/f'inputs/pairs_{name}.txt').open('w') as f:
 f.write(str(len(out))+'\n')
 for d in out:
  fs=[d['families'][i] for i in d['keep']];f.write(f"{d['id']} 1 {len(fs)}\n")
  for x in fs:
   f.write(str(x['free'])+'\n'+' '.join(str(v) for e in x['row'] for v in e)+'\n'+str(len(x['eq']))+'\n')
   for e in x['eq']:f.write(' '.join(map(str,e))+'\n')
