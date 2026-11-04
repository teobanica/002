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
sys.path.insert(0,str(R.parent/'proofs/AAA_classes11_19/code'))
import torus3 as T

class Lat:
 def __init__(self,eq=(),mod=120):
  self.mod=mod;self.h=[[0]*4 for _ in range(4)];self.h[3][3]=mod;self.ok=True
  for e in eq:
   if not self.insert(e):break
 def insert(self,v):
  v=list(v)
  for j in range(4):
   while v[j]:
    if not self.h[j][j]:
     if v[j]<0:v=[-x for x in v]
     self.h[j]=v;v=[0]*4;break
    q=v[j]//self.h[j][j];v=[x-q*y for x,y in zip(v,self.h[j])]
    if v[j]:v,self.h[j]=self.h[j],v
  if self.h[3][3]!=self.mod:self.ok=False;return False
  for j in range(4):
   if self.h[j][j]:
    for i in range(j):
     if self.h[i][i]:
      q=self.h[i][j]//self.h[j][j];self.h[i]=[x-q*y for x,y in zip(self.h[i],self.h[j])]
  return True
 def implies(self,v):
  v=list(v)
  for j in range(4):
   if self.h[j][j]:
    q=v[j]//self.h[j][j];v=[x-q*y for x,y in zip(v,self.h[j])]
  return not any(v)
 def eq(self):return [e for e in self.h[:3] if any(e)]
 def key(self):return tuple(map(tuple,self.eq()))
 def subset(self,b):return all(self.implies(e) for e in b.eq())

def covers(a,b,L):
 if a['free']>b['free']:return False
 values=[];w=4+a['free']
 for k in range(b['free']):
  j=next(j for j,e in enumerate(b['row']) if e[4+k])
  assert b['row'][j][4+k]==1 and sum(b['row'][j][4:])==1
  values.append([a['row'][j][u]-b['row'][j][u] for u in range(4)]+a['row'][j][4:])
 for x,y in zip(a['row'],b['row']):
  v=[y[u]-x[u] for u in range(4)]+[-z for z in x[4:]]
  for k,z in enumerate(values):
   for u in range(w):v[u]+=y[4+k]*z[u]
  if any(v[4:]) or not L.implies(v[1:4]+[v[0]*(L.mod//30)]):return False
 return True

def specialize(row,t):
 return [[(F(e[0],30)+sum(F(e[i+1])*t['o'][i] for i in range(3)))%1]+[sum(e[i+1]*t['V'][i][j] for i in range(3)) for j in range(t['dim'])]+e[4:] for e in row]

data=json.loads((R/'outputs/full_rows_11_19.json').read_text())
seeds={s['class_id']:s for s in json.loads((R.parent/'proofs/AAA_classes11_19/inputs/components.json').read_text())}
roots=[]
for name in ['good_planes','good_lines']:
 roots+=json.loads((R.parent/f'proofs/AAA_classes11_19/outputs/{name}.json').read_text())
strata=[];edges=[];done={};stats=[]
for d in data:
 cid=d['id'];fs=d['families'];groups=defaultdict(list)
 for i,f in enumerate(fs):groups[tuple(map(tuple,f['eq']))].append(i)
 def visit(t):
  key=(cid,T.key(t))
  if key in done:return done[key]
  sid=len(strata);done[key]=sid;strata.append(None)
  M=lcm(120,*(c.denominator for c in t['c']));L=Lat(T.equations(t,M),M);assert L.ok
  generic=[];proper=[]
  for eq,indices in groups.items():
   scaled=[list(e[:-1])+[e[-1]*(M//30)] for e in eq]
   if all(L.implies(e) for e in scaled):generic+=indices
   else:
    joint=Lat(L.eq()+scaled,M)
    if joint.ok:proper.append((joint,indices))
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
  strata[sid]=dict(id=sid,class_id=cid,torus=T.enc(t),base=t['dim'],rows=[specialize(row,t) for row in seeds[cid]['rows']],families=families)
  print('stratum',sid,'class',cid,'dim',t['dim'],'families',len(keep),'child loci',len(maximal),flush=True)
  for child in maximal:
   for u in T.split(child.eq(),M):
    assert u['dim']<t['dim'] and T.subset(u,t)
    target=visit(u);edges.append([sid,target])
  return sid
 rootids=[visit(T.decode(x['torus'])) for x in roots if x['class_id']==cid]
 stats.append(dict(class_id=cid,root_strata=rootids))
# Fractions are stored as numerator/denominator pairs, integer coefficients stay integers.
(R/'outputs/strata_11_19.json').write_text(json.dumps(T.enc(dict(strata=strata,edges=edges,roots=stats)),separators=(',',':'))+'\n')
print('Finished',len(strata),'strata',dict(Counter(s['base'] for s in strata)),flush=True)
