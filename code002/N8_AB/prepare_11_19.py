from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
from math import lcm
from collections import Counter
import json
R=Path(__file__).resolve().parents[1]
W=R/'stage11_19'
for sub in ['inputs','outputs','checks']:(W/sub).mkdir(parents=True,exist_ok=True)
def canonical(H):
 best=None;cert=None
 for rp in permutations(range(3)):
  for p in range(8):
   cols=[tuple((H[r][j]-H[r][p]-H[rp[0]][j]+H[rp[0]][p])%1 for r in rp) for j in range(8)]
   cp=sorted(range(8),key=lambda j:cols[j]);key=tuple(cols[j] for j in cp)
   if best is None or key<best:best=key;cert=dict(rows=rp,pivot_column=p,columns=cp)
 return best,cert
old=json.loads((R/'outputs/strata.json').read_text());point_reps={}
for s in old['strata']:
 if s['base']==0:
  key,cert=canonical([[F(e[0],old['modulus']) for e in row] for row in s['rows']]);point_reps[key]=dict(stage='21_22',id=s['id'],transform=cert)
data=json.loads((R/'outputs/strata_11_19.json').read_text());keep=[];maps=[];mod=60
for s in data['strata']:
 if s['base']==0:
  key,cert=canonical([[F(*e[0]) for e in row] for row in s['rows']])
  if key in point_reps:
   maps.append(dict(source=s['id'],target=point_reps[key],transform=cert));continue
  point_reps[key]=dict(stage='11_19',id=len(keep),transform=cert)
 k=dict(s,id=len(keep),raw_stratum=s['id']);keep.append(k)
 for row in k['rows']+[f['row'] for f in k['families']]:
  for e in row:mod=lcm(mod,F(*e[0]).denominator)
for s in keep:
 for row in s['rows']+[f['row'] for f in s['families']]:
  for e in row:e[0]=int(F(*e[0])*mod)
(W/'outputs/strata.json').write_text(json.dumps(dict(modulus=mod,strata=keep,point_equivalences=maps),separators=(',',':'))+'\n')
with (W/'inputs/pairs_strata.txt').open('w') as f:
 f.write(str(len(keep))+'\n')
 for s in keep:
  f.write(f"{s['id']} {s['base']} {len(s['families'])}\n")
  for x in s['families']:f.write(str(x['free'])+'\n'+' '.join(str(v) for e in x['row'] for v in e)+'\n0\n')
print('Retained',len(keep),'strata',dict(Counter(s['base'] for s in keep)),'modulus',mod,'point equivalences',len(maps))
