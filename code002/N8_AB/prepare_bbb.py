"""All normalized label lifts and Smith components of the 39 BBB support classes."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,permutations
from math import lcm
import json,sys
import sympy as S
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
R=Path(__file__).resolve().parents[1];W=R/'stage_bbb'
for name in ['inputs','outputs','checks']:(W/name).mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from exact_rows import profiles,quotient
cs=json.loads((R.parent/'proofs/consolidation/outputs/canonical.json').read_text())['triples']['BBB']['classes']
def labels(part):
 out=[]
 for choices in product(*(list(permutations(range(1,len(b)))) for b in part)):
  row=[0]*8
  for b,pr in zip(part,choices):
   for j,k in enumerate(pr,1):row[b[j]]=30*k//len(b)
  out.append(row)
 return out
out=[];counts=[];mod=30;skipped=0
for c in cs:
 parts=[[[0,1,2],[3,4,5],[6,7]],c['Q'],c['R']];A=S.zeros(8,7)
 for k,part in enumerate(parts):
  for i,b in enumerate(part):
   if k<2 and i==0:continue
   ix=[0,2,4][k]+i-(k<2)
   for j in b:A[j,ix]=[1,-1,-1][k]
 sd,ud,vd=smith_normal_decomp(DomainMatrix.from_Matrix(A).convert_to(ZZ));D,U,V=[x.to_Matrix() for x in (sd,ud,vd)];rank=A.rank();dim=7-rank
 assert dim==c['dim'] and U*A*V==D
 lifts=branches=0
 for labs in product(*(labels(p) for p in parts)):
  rhs=S.Matrix([S.Rational(-labs[0][j]+labs[1][j]+labs[2][j],30) for j in range(8)]);b=U*rhs
  if any(x.q!=1 for x in b[rank:]):continue
  lifts+=1
  for ks in product(*(range(abs(int(D[i,i]))) for i in range(rank))):
   y=S.zeros(7,1)
   for i in range(rank):y[i]=(b[i]+ks[i])/D[i,i]
   offset=V*y;basis=V[:,rank:];rows=[[[F(0)]+[0]*dim for _ in range(8)]]
   for k in range(2):
    row=[]
    for j in range(8):
     bi=next(i for i,bb in enumerate(parts[k]) if j in bb);ix=2*k+bi-1
     v=[F(labs[k][j],30)+(F(offset[ix]) if bi else 0)]+([int(z) for z in basis[ix,:]] if bi else [0]*dim);row.append(v)
    anchor=row[0];rows.append([[(e[0]-anchor[0])%1]+[x-y for x,y in zip(e[1:],anchor[1:])] for e in row])
   ps=[profiles(quotient(rows[i],rows[j])) for i,j in [(0,1),(0,2),(1,2)]];assert all('B' in p for p in ps)
   branches+=1
   if all('A' in p for p in ps):skipped+=1;continue
   for row in rows:
    for e in row:mod=lcm(mod,e[0].denominator)
   out.append(dict(class_id=c['id'],label_lift=lifts-1,torsion=list(ks),dim=dim,rows=rows,profiles=ps))
 assert lifts==c['normalized_label_lifts'],(c['id'],lifts,c['normalized_label_lifts'])
 counts.append(dict(class_id=c['id'],label_lifts=lifts,smith_branches=branches))
 print(c['id'],dim,lifts,branches,flush=True)
for i,s in enumerate(out):
 s['id']=i
 for row in s['rows']:
  for e in row:e[0]=int(e[0]*mod)
(W/'inputs/components.json').write_text(json.dumps(dict(modulus=mod,components=out,counts=counts,uniformly_AAA_skipped=skipped),separators=(',',':'))+'\n')
print('All components',sum(c['smith_branches'] for c in counts),'retained',len(out),'uniformly AAA',skipped,'modulus',mod,flush=True)
