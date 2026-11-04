"""Independent integer verification of every reported fourth-row support witness."""
from pathlib import Path
from itertools import combinations,permutations,product
import json
R=Path(__file__).resolve().parents[1]
def matchings(left):
 if not left:yield [];return
 a=left[0]
 for b in left[1:]:
  for rest in matchings([x for x in left if x not in (a,b)]):yield [(a,b)]+rest
supports=[]
for m in matchings(list(range(8))):supports.append([(a,b,15) for a,b in m])
for a,b in combinations(range(8),2):
 r=[j for j in range(8) if j not in (a,b)]
 for mask in range(64):
  if mask.bit_count()!=3:continue
  x=[r[j] for j in range(6) if mask>>j&1];y=[r[j] for j in range(6) if not(mask>>j&1)]
  if x>=y:continue
  for p,q in product(permutations([1,2]),repeat=2):supports.append([(x[0],x[j],10*p[j-1]) for j in (1,2)]+[(y[0],y[j],10*q[j-1]) for j in (1,2)]+[(a,b,15)])
assert len(supports)==1225

def implies(eq,v,mod=30):
 # Input rows have coefficients first and mod-30 constant last.
 for e in eq:
  j=next(i for i,x in enumerate(e[:-1]) if x)
  q=v[j]//e[j];v=[a-q*b for a,b in zip(v,e)]
 return not any(v[:-1]) and v[-1]%mod==0

def verify(data,seeds,mod=30):
 count=0
 for d,s in zip(data,seeds):
  assert d['id']==s['id'];base=d['dim']
  for f in d['families']:
   for r,sid in zip(s['rows'],f['supports']):
    for a,b,c in supports[sid]:
     va,vb=f['row'][a],f['row'][b];e=[vb[k]-va[k]-(r[b][k]-r[a][k] if k<=base else 0) for k in range(len(va))];e[0]-=c*(mod//30)
     assert not any(e[base+1:])
     assert implies(f['eq'],e[1:base+1]+[e[0]],mod),(d['id'],f['supports'],a,b,e)
   count+=1
 return count
if __name__=='__main__':
 seeds=json.loads((R/'inputs/circle_components.json').read_text());data=json.loads((R/'outputs/full_rows_all.json').read_text())
 n=verify(data,seeds)
 summary=dict(circle_row_witnesses=n,all_circle_witnesses_valid=True)
 if (R/'outputs/full_rows_11_19.json').exists():
  ds=json.loads((R/'outputs/full_rows_11_19.json').read_text());ss=json.loads((R.parent/'proofs/AAA_classes11_19/inputs/components.json').read_text())
  ss=[dict(s,id=s['class_id']) for s in ss];summary['classes11_19_row_witnesses']=verify(ds,ss)
 (R/'checks/row_witness_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(summary)
