"""Exact Smith certificates, including unimodularity, for every leaf matrix."""
from pathlib import Path
import json,sys
import sympy as S
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
W=Path(sys.argv[1]).resolve();ss=json.loads((W/'outputs/strata.json').read_text())['strata'];keys=set()
for s in ss:
 base=s['base'];fs=s['families'];assert S.Matrix([e[1:] for row in s['rows'] for e in row]).rank()==base
 for f in fs:
  for k in range(f['free']):assert any(e[1+base+k]==1 and sum(abs(v) for v in e[1+base:])==1 for e in f['row'])
 raw=json.loads((W/f"outputs/cliques_stratum_{s['id']:02d}.json").read_text())[0]['squares']
 for loc in raw:
  n=base+sum(fs[i]['free'] for i in loc['families']);keys.add((n,tuple(tuple(e[:-1]) for e in loc['eq'])))
certs=[]
for n,key in sorted(keys):
 if not key:certs.append(dict(variables=n,A=[],D=[],U=[],V=S.eye(n).tolist()));continue
 A=S.Matrix(key);sd,ud,vd=smith_normal_decomp(DomainMatrix.from_Matrix(A).convert_to(ZZ));D,U,V=[x.to_Matrix() for x in (sd,ud,vd)]
 assert U*A*V==D and abs(U.det())==abs(V.det())==1
 for i in range(D.rows):
  for j in range(D.cols):assert i==j or D[i,j]==0
 ds=[abs(int(D[i,i])) for i in range(min(D.rows,D.cols)) if D[i,i]];assert all(b%a==0 for a,b in zip(ds,ds[1:]))
 certs.append(dict(variables=n,A=[list(r) for r in key],D=[[int(x) for x in r] for r in D.tolist()],U=[[int(x) for x in r] for r in U.tolist()],V=[[int(x) for x in r] for r in V.tolist()]))
p=W/'outputs/smith_certificates.json';tmp=p.with_suffix('.pending');tmp.write_text(json.dumps(certs,separators=(',',':'))+'\n');tmp.replace(p)
summary=dict(all_checks_passed=True,coefficient_matrices=len(certs),all_U_V_unimodular=True,all_stratum_parameter_maps_have_full_rank=True,all_row_phase_coordinates_independent=True)
(W/'checks/smith_verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
