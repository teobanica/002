from pathlib import Path
from fractions import Fraction as F
from math import lcm
from collections import Counter
import json,sys
W=Path(sys.argv[1]).resolve();source=sys.argv[2]
doc=json.loads((W/'outputs'/source).read_text());out=[];mod=60
for s in doc['strata']:
 if 'families' not in s:continue
 t=dict(s,id=len(out),raw_stratum=s['id']);out.append(t)
 for row in t['rows']+[f['row'] for f in t['families']]:
  for e in row:mod=lcm(mod,F(*e[0]).denominator)
for s in out:
 for row in s['rows']+[f['row'] for f in s['families']]:
  for e in row:e[0]=int(F(*e[0])*mod)
 assert max((a['free']+b['free']+s['base'] for a in s['families'] for b in s['families']),default=0)<9
(W/'outputs/strata.json').write_text(json.dumps(dict(modulus=mod,strata=out,raw_to_search={s['raw_stratum']:s['id'] for s in out}),separators=(',',':'))+'\n')
with (W/'inputs/pairs_strata.txt').open('w') as f:
 f.write(str(len(out))+'\n')
 for s in out:
  f.write(f"{s['id']} {s['base']} {len(s['families'])}\n")
  for x in s['families']:f.write(str(x['free'])+'\n'+' '.join(str(v) for e in x['row'] for v in e)+'\n0\n')
print(len(out),'explicit strata',dict(Counter(s['base'] for s in out)),'modulus',mod)
