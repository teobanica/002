"""Cached Smith extraction using integer affine images, with exact cycle checks.

Caches depend only on integer coefficient matrices. Every torsion branch is
retained. Duplicate leaf systems are recorded explicitly before elimination.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import product,combinations
from collections import Counter
from math import gcd,lcm
import json,sys
import sympy as S
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
HOME_ROOT=Path(__file__).resolve().parents[1];W=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(HOME_ROOT.parent/'n8_aaa_complete/code'))
from exact_rows import profiles,quotient
from prepare_affine import write_input
D=json.loads((W/'outputs/strata.json').read_text());MOD=D['modulus'];strata=D['strata']

@lru_cache(maxsize=20000)
def smith(Akey,n):
 m=len(Akey)
 if not m:return n,1,tuple(tuple(int(i==j) for j in range(n)) for i in range(n)),tuple(() for _ in range(n)),(tuple(0 for _ in range(n)),),()
 A=S.Matrix(Akey);sd,ud,vd=smith_normal_decomp(DomainMatrix.from_Matrix(A).convert_to(ZZ));diag,U,V=[x.to_Matrix() for x in (sd,ud,vd)];assert U*A*V==diag
 rank=sum(diag[i,i]!=0 for i in range(min(m,n)));div=lcm(1,*(abs(int(diag[i,i])) for i in range(rank)))
 u=[[int(x) for x in row] for row in U.tolist()];v=[[int(x) for x in row] for row in V.tolist()]
 weights=tuple(tuple(-sum(v[i][k]*(div//int(diag[k,k]))*u[k][j] for k in range(rank)) for j in range(m)) for i in range(n))
 torsion=[]
 for ks in product(*(range(abs(int(diag[k,k]))) for k in range(rank))):torsion.append(tuple(sum(v[i][k]*(div//int(diag[k,k]))*ks[k] for k in range(rank)) for i in range(n)))
 return n-rank,div,tuple(tuple(row[rank:]) for row in v),weights,tuple(torsion),tuple(tuple(row) for row in u[rank:])

@lru_cache(maxsize=200000)
def check(q):
 p=profiles(q);assert 'A' in p or 'B' in p,(q,p)
 return p

def decoded(key):return tuple(tuple((F(e[0],e[1]),)+tuple(map(F,e[2:])) for e in row) for row in key)
def encoded(H):return [[[[x.numerator,x.denominator] for x in e] for e in row] for row in H]
seen={};out=[];counts=[];duplicates=[];modulus=MOD
for s in strata:
 sid=s['id'];raw=json.loads((W/f'outputs/cliques_stratum_{sid:02d}.json').read_text())[0]['squares'];fs=s['families'];base=s['base'];unique={};raw_branches=unique_branches=0
 for li,loc in enumerate(raw):
  ids=tuple(loc['families']);eq=tuple(map(tuple,loc['eq']));chosen=[fs[i] for i in ids];nf=sum(f['free'] for f in chosen);n=nf+base
  Akey=tuple(e[:-1] for e in eq);dim,div,basis,weights,torsion,redundant=smith(Akey,n);assert dim<=5,(sid,li,dim)
  raw_branches+=len(torsion);lk=(ids,eq)
  if lk in unique:duplicates.append([sid,li,unique[lk]]);continue
  unique[lk]=li;unique_branches+=len(torsion)
  const=[e[-1] for e in eq]
  assert all(sum(a*b for a,b in zip(row,const))%MOD==0 for row in redundant)
  offset=[sum(a*b for a,b in zip(row,const)) for row in weights];den=MOD*div
  # Sparse maps from all latent variables to each of the 64 matrix entries.
  maps=[]
  for row in s['rows']:
   maps.append([(e[0],[(nf+k,e[1+k]) for k in range(base) if e[1+k]]) for e in row])
  pos=0
  for f in chosen:
   maps.append([(e[0],[(nf+k,e[1+k]) for k in range(base) if e[1+k]]+[(pos+k,e[1+base+k]) for k in range(f['free']) if e[1+base+k]]) for e in f['row']]);pos+=f['free']
  affine=[]
  for row in maps:
   affine.append([(c*div+sum(k*offset[i] for i,k in sparse),tuple(sum(k*basis[i][a] for i,k in sparse) for a in range(dim)),sparse) for c,sparse in row])
  for bi,t in enumerate(torsion):
   H=[]
   for row in affine:
    rr=[]
    for c,slopes,sparse in row:
     num=(c+MOD*sum(k*t[i] for i,k in sparse))%den;g=gcd(num,den);rr.append((num//g,den//g)+slopes)
    H.append(tuple(rr))
   assert all(e[0]==0 and not any(e[2:]) for e in H[0]);assert all(row[0][0]==0 and not any(row[0][2:]) for row in H)
   key=tuple(H[:3])+tuple(sorted(H[3:]));origin=dict(stratum=sid,locus=li,torsion_branch=bi)
   if key in seen:out[seen[key]]['origins'].append(origin);continue
   hh=decoded(key);pp=[check(tuple(tuple(e) for e in quotient(hh[i],hh[j]))) for i,j in combinations(range(8),2)]
   for row in key:
    for e in row:modulus=lcm(modulus,e[1])
   seen[key]=len(out);out.append(dict(square=encoded(hh),dimension=dim,origins=[origin],profiles=pp))
 counts.append(dict(stratum=sid,loci=len(raw),distinct_leaf_systems=len(unique),smith_branches=raw_branches,distinct_leaf_smith_branches=unique_branches))
 print(sid,'loci',len(raw),'distinct',len(unique),'branches',raw_branches,'unique squares',len(out),flush=True)
# Validate complete serialization before publishing outputs.
for name,data in [('completion_squares.json',out),('duplicate_leaf_systems.json',duplicates),('square_extraction_summary.json',dict(modulus=modulus,counts=counts,unique=len(out),dimensions=dict(Counter(f['dimension'] for f in out)),all_strata_finished=True,smith_cache=str(smith.cache_info())))]:
 p=W/'outputs'/name;blob=json.dumps(data,separators=(',',':'))+'\n';json.loads(blob);tmp=p.with_suffix('.pending');tmp.write_text(blob);tmp.replace(p)
targets=json.loads((HOME_ROOT.parent/'n8_aaa_complete/outputs/ten_representatives.json').read_text());write_input(W/'inputs/completion_embeddings.txt',targets,out,modulus)
print('Finished',len(out),'families, modulus',modulus,dict(Counter(f['dimension'] for f in out)),smith.cache_info(),flush=True)
