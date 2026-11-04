"""Check saved triple equivalences directly, without repeating their searches."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json
R=Path(__file__).resolve().parents[1];W=R/'stage_uniform'
old={}
for stage,path in [('21_22',R),('11_19',R/'stage11_19')]:
 d=json.loads((path/'outputs/strata.json').read_text());M=d['modulus'];old[stage]={s['id']:[[(F(e[0],M),)+tuple(e[1:]) for e in row] for row in s['rows']] for s in d['strata']}
doc=json.loads((W/'outputs/strata_equivalent.json').read_text());old['uniform']={s['id']:[[(F(*e[0]),)+tuple(e[1:]) for e in row] for row in s['rows']] for s in doc['strata']}
def image(H,w):
 p=w['rows'];q=w['columns'];a=w['pivot_column'];assert sorted(p)==list(range(3)) and sorted(q)==list(range(8)) and 0<=a<8
 dim=len(H[0][0])-1;out=[]
 for j in q:
  col=[]
  for r in p:
   c=(H[r][j][0]-H[r][a][0]-H[p[0]][j][0]+H[p[0]][a][0])%1
   if dim:
    assert dim==1;g=w['slope_gcd'];sgn=w['parameter_sign'];shift=F(*w['parameter_shift']);assert g>0 and sgn in [-1,1]
    k=H[r][j][1]-H[r][a][1]-H[p[0]][j][1]+H[p[0]][a][1];assert k%g==0;k=k//g*sgn;col.append(((c+k*shift)%1,k))
   else:col.append((c,))
  out.append(tuple(col))
 return tuple(out)
counts=Counter()
for s in doc['strata']:
 if 'alias' not in s:continue
 a=s['alias'];source=old['uniform'][s['id']];target=old[a['stage']][a['id']]
 assert image(source,s['transform'])==image(target,a['transform']),(s['id'],a)
 if a['stage']=='uniform':assert a['id']<s['id']
 counts[(s['base'],a['stage'])]+=1
summary=dict(all_equivalences_verified=True,total=sum(counts.values()),counts={str(k):v for k,v in counts.items()})
(W/'checks/alias_verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
