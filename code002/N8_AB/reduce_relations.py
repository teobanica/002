"""Remove pair loci contained in other pair loci, with exact lattice witnesses."""
from pathlib import Path
import sys,json
HOME_ROOT=Path(__file__).resolve().parents[1]
R=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HOME_ROOT
data=json.loads((R/'outputs/pairs_strata.json').read_text())
MOD=json.loads((R/'outputs/strata.json').read_text())['modulus']
def implies(eq,v):
 v=v.copy()
 for e in eq:
  j=next(i for i,x in enumerate(e[:-1]) if x)
  q=v[j]//e[j]
  v=[a-q*b for a,b in zip(v,e)]
 return not any(v[:-1]) and v[-1]%MOD==0
def subset(a,b):return all(implies(a,e) for e in b)
for d in data:
 before=after=0
 for p in d['pairs']:
  rs=p['relations'];keep=[];maps=[];before+=len(rs)
  for i in sorted(range(len(rs)),key=lambda i:len(rs[i]['eq'])):
   for j in keep:
    if subset(rs[i]['eq'],rs[j]['eq']):maps.append([i,j]);break
   else:keep.append(i)
  # Equal-rank sublattice inclusions require a second pass in either direction.
  final=[]
  for i in keep:
   if any(j!=i and subset(rs[i]['eq'],rs[j]['eq']) and (not subset(rs[j]['eq'],rs[i]['eq']) or j<i) for j in keep):continue
   final.append(i)
  for i in range(len(rs)):
   assert any(subset(rs[i]['eq'],rs[j]['eq']) for j in final)
  p['relations']=[dict(rs[i],raw_index=i) for i in final];p['raw_relation_count']=len(rs);after+=len(final)
 print(d['class_id'],before,'->',after,flush=True)
(R/'outputs/pairs_strata_maximal.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
