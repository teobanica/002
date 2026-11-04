"""Complete base stratification for all completions on the known good loci.

At each connected base torus, retain every generic fourth-row family. A proper
base locus needs further treatment only if it carries a row outside this list.
Maximal such loci are Smith-split and treated recursively. All computations are
integer lattice operations and exact Fraction substitutions.
"""
from pathlib import Path
import sys,json
from fractions import Fraction as F
from math import lcm
from collections import defaultdict,Counter
R=Path(__file__).resolve().parents[1]
import torus_nd as T
from itertools import permutations
N=3

class Lat:
 def __init__(self,eq=(),mod=120):
  self.mod=mod;self.h=[[0]*(N+1) for _ in range((N+1))];self.h[N][N]=mod;self.ok=True
  for e in eq:
   if not self.insert(e):break
 def insert(self,v):
  v=list(v)
  for j in range((N+1)):
   while v[j]:
    if not self.h[j][j]:
     if v[j]<0:v=[-x for x in v]
     self.h[j]=v;v=[0]*(N+1);break
    q=v[j]//self.h[j][j];v=[x-q*y for x,y in zip(v,self.h[j])]
    if v[j]:v,self.h[j]=self.h[j],v
  if self.h[N][N]!=self.mod:self.ok=False;return False
  for j in range((N+1)):
   if self.h[j][j]:
    for i in range(j):
     if self.h[i][i]:
      q=self.h[i][j]//self.h[j][j];self.h[i]=[x-q*y for x,y in zip(self.h[i],self.h[j])]
  return True
 def implies(self,v):
  v=list(v)
  for j in range((N+1)):
   if self.h[j][j]:
    q=v[j]//self.h[j][j];v=[x-q*y for x,y in zip(v,self.h[j])]
  return not any(v)
 def eq(self):return [e for e in self.h[:N] if any(e)]
 def key(self):return tuple(map(tuple,self.eq()))
 def subset(self,b):return all(self.implies(e) for e in b.eq())

def covers(a,b,L):
 if a['free']>b['free']:return False
 values=[];w=(N+1)+a['free']
 for k in range(b['free']):
  j=next(j for j,e in enumerate(b['row']) if e[(N+1)+k])
  assert b['row'][j][(N+1)+k]==1 and sum(b['row'][j][(N+1):])==1
  values.append([a['row'][j][u]-b['row'][j][u] for u in range((N+1))]+a['row'][j][(N+1):])
 for x,y in zip(a['row'],b['row']):
  v=[y[u]-x[u] for u in range((N+1))]+[-z for z in x[(N+1):]]
  for k,z in enumerate(values):
   for u in range(w):v[u]+=y[(N+1)+k]*z[u]
  if any(v[(N+1):]) or not L.implies(v[1:(N+1)]+[v[0]*(L.mod//P)]):return False
 return True

def specialize(row,t):
 return [[(F(e[0],P)+sum(F(e[i+1])*t['o'][i] for i in range(N)))%1]+[sum(e[i+1]*t['V'][i][j] for i in range(N)) for j in range(t['dim'])]+e[(N+1):] for e in row]


W=R/'stage_bbb'
document=json.loads((W/'inputs/reduced_components.json').read_text());P=document['modulus'];components=document['components']
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from exact_rows import profiles,quotient
point_reps={}
def point_key(rows):
    M=lcm(1,*(e[0].denominator for row in rows for e in row))
    H=[[int(e[0]*M) for e in row] for row in rows];best=None;cert=None
    for rp in permutations(range(3)):
        for p in range(8):
            cols=[tuple((H[r][j]-H[r][p]-H[rp[0]][j]+H[rp[0]][p])%M for r in rp) for j in range(8)]
            cp=sorted(range(8),key=lambda j:cols[j]);key=tuple(cols[j] for j in cp)
            if best is None or key<best:best=key;cert=dict(rows=rp,pivot_column=p,columns=cp)
    return (M,best),cert
for stage,path in [('21_22',R),('11_19',R/'stage11_19')]:
    doc=json.loads((path/'outputs/strata.json').read_text());M=doc['modulus']
    for s in doc['strata']:
        if s['base']!=0:continue
        key,cert=point_key([[[F(e[0],M)] for e in row] for row in s['rows']])
        point_reps[key]=dict(stage=stage,id=s['id'],transform=cert)
strata=[];edges=[];done={};stats=[]
for seed in components:
 N=T.N=seed['dim'];cid=seed['component_index']
 d=json.loads((W/f'outputs/full_rows_{cid:02d}.json').read_text())[0]
 fs=d['families'];groups=defaultdict(list)
 for i,f in enumerate(fs):groups[tuple(map(tuple,f['eq']))].append(i)
 def visit(t,active=None):
  key=(cid,T.key(t))
  if key in done:return done[key]
  sid=len(strata);done[key]=sid;strata.append(None)
  rows=[specialize(row,t) for row in seed['rows']]
  if any('A' in profiles(quotient(rows[i],rows[j])) for i,j in [(0,1),(0,2),(1,2)]):
   strata[sid]=dict(id=sid,component=cid,class_id=seed['class_id'],base=t['dim'],torus=T.enc(t),rows=rows,excluded='base A edge')
   return sid
  if t['dim']==0:
   pk,cert=point_key(rows)
   if pk in point_reps:
    strata[sid]=dict(id=sid,component=cid,class_id=seed['class_id'],base=0,torus=T.enc(t),rows=rows,alias=point_reps[pk],transform=cert)
    return sid
   point_reps[pk]=dict(stage='uniform',id=sid,transform=cert)
  M=lcm(120,*(c.denominator for c in t['c']));L=Lat(T.equations(t,M),M);assert L.ok
  generic=[];proper=[];alive={}
  if active is None:active=groups
  for eq,indices in active.items():
   scaled=[list(e[:-1])+[e[-1]*(M//P)] for e in eq]
   if all(L.implies(e) for e in scaled):generic+=indices;alive[eq]=indices
   else:
    joint=Lat(L.eq()+scaled,M)
    if joint.ok:proper.append((joint,indices));alive[eq]=indices
  keep=[]
  for i in sorted(generic,key=lambda i:(-fs[i]['free'],i)):
   if not any(covers(fs[i],fs[j],L) for j in keep):keep.append(i)
  children={}
  for joint,indices in proper:
   if any(not any(covers(fs[i],fs[j],joint) for j in keep) for i in indices):children[joint.key()]=joint
  maximal=[]
  for child in sorted(children.values(),key=lambda x:len(x.eq())):
   if any(child.subset(other) for other in maximal):continue
   maximal=[other for other in maximal if not other.subset(child)];maximal.append(child)
  families=[dict(raw_index=i,free=fs[i]['free'],row=specialize(fs[i]['row'],t)) for i in keep]
  strata[sid]=dict(id=sid,component=cid,class_id=seed['class_id'],torus=T.enc(t),base=t['dim'],rows=rows,families=families)
  if t['dim']>1 or sid%100==0:print('stratum',sid,'component',cid,'dim',t['dim'],'families',len(keep),'child loci',len(maximal),flush=True)
  for child in maximal:
   for u in T.split(child.eq(),M):
    assert u['dim']<t['dim'] and T.subset(u,t)
    target=visit(u,alive);edges.append([sid,target])
  return sid
 root=visit(T.torus([],[]))
 stats.append(dict(component=cid,root=root))
 snapshot=T.enc(dict(strata=strata,edges=edges,roots=stats,completed_components=len(stats)))
 pending=W/'outputs/strata_checkpoint.pending';pending.write_text(json.dumps(snapshot,separators=(',',':'))+'\n');pending.replace(W/'outputs/strata_checkpoint.json')
 print('component',cid,'finished, cumulative',len(strata),'strata, explicit',sum('alias' not in s for s in strata),flush=True)
(W/'outputs/strata_bbb.json').write_text(json.dumps(T.enc(dict(strata=strata,edges=edges,roots=stats)),separators=(',',':'))+'\n')
print('Finished',len(strata),'strata',dict(Counter(s['base'] for s in strata)),flush=True)
