"""Independent cycle generation, graph traversal and rational-minor solving.
Compares complete maximal joint solution loci with the primary census.
"""
from itertools import permutations
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from collections import Counter
import json,time,sys
from families import family
from minor_solver import solve as minor_solve
from census import canonical,contains,deserialize
ROOT=Path(__file__).resolve().parent

def patterns():
    out=set()
    for p in permutations(range(7)):
        edges=[]
        for b in (p[:3],p[3:5],p[5:]):
            r=min(b);a=b.index(r);m=len(b)
            edges.extend((j,r,((k-a)%m)*(6//m)) for k,j in enumerate(b) if j!=r)
        out.add(tuple(sorted(edges)))
    return sorted(out)
PAT=patterns();assert len(PAT)==210

def egcd(a,b):
    if b==0:return abs(a),1 if a>=0 else -1,0
    g,x,y=egcd(b,a%b);return g,y,x-(a//b)*y

def kernel(A,n):
    K=[[int(i==j) for j in range(n)] for i in range(n)];d=n
    for a in A:
        r=[sum(a[i]*K[i][j] for i in range(n)) for j in range(d)]
        if not any(r):continue
        for j in range(1,d):
            if r[j]==0:continue
            x,y=r[0],r[j];g,s,t=egcd(x,y)
            for row in K:
                u,v=row[0],row[j];row[0]=s*u+t*v;row[j]=-(y//g)*u+(x//g)*v
            r[0]=g;r[j]=0
        K=[row[1:] for row in K];d-=1
    return tuple(tuple(r) for r in K)

@lru_cache(None)
def saturated(B):
    n=len(B);L=kernel(tuple(zip(*B)),n)
    return kernel(tuple(zip(*L)),n)

def canon(c,B):
    B=tuple(tuple(v) for v in B)
    return canonical(c,saturated(B))

@lru_cache(None)
def solve(A,b,n):
    if n==0:return (((),()),) if all(v%1==0 for v in b) else ()
    return tuple(minor_solve(A,b,n))

def edges(p,c,B):return [(i,j,(e+c[i]-c[j])%6,tuple(x-y for x,y in zip(B[i],B[j]))) for i,j,e in p]

def joint(eq,d):
    adj=[[] for _ in range(7)]
    for i,j,c,B in eq:
        adj[j].append((i,c,B));adj[i].append((j,-c,tuple(-v for v in B)))
    labs=[None]*7;conditions={}
    for root in range(7):
        if labs[root] is not None:continue
        labs[root]=(root,0,(0,)*d);todo=[root]
        while todo:
            j=todo.pop();_,c0,B0=labs[j]
            for i,c,B in adj[j]:
                c=(c+c0)%6;B=tuple(x+y for x,y in zip(B,B0))
                if labs[i] is None:labs[i]=(root,c,B);todo.append(i)
                else:
                    _,ci,Bi=labs[i];dc=(c-ci)%6;db=tuple(x-y for x,y in zip(B,Bi))
                    if not any(db):
                        if dc:return
                    else:
                        if next(v for v in db if v)<0:dc=(-dc)%6;db=tuple(-v for v in db)
                        if db in conditions and conditions[db]!=dc:return
                        conditions[db]=dc
    cs=sorted(conditions.items());A=tuple(b for b,c in cs);rhs=tuple(Q(-c,6)%1 for b,c in cs)
    free=sorted({r for r,c,b in labs}-{0})
    for c,B in solve(A,rhs,d):
        k=len(B[0]);cc=list(c);BB=[row+(0,)*len(free) for row in B]
        for r,v,b in labs:
            cc.append((Q(v,6)+sum(x*y for x,y in zip(b,c)))%1)
            BB.append(tuple(sum(b[i]*B[i][j] for i in range(d)) for j in range(k))+tuple(int(r==s) for s in free))
        yield canon(cc,BB)

def restrict(comp,pat,s):
    c,B,piv=comp;d=s['parameters'];A=[];rhs=[]
    for i,j,e in pat:
        coeff=[0]*(d+7);coeff[d+i]=1;coeff[d+j]=-1
        for h in range(d):coeff[h]-=s['basis'][2][i][h]-s['basis'][2][j][h]
        ar=tuple(sum(coeff[h]*B[h][k] for h in range(d+7)) for k in range(len(piv)))
        br=(Q(e+s['constants'][2][i]-s['constants'][2][j],6)-sum(coeff[h]*c[h] for h in range(d+7)))%1
        if not any(ar):
            if br:return
        else:A.append(ar);rhs.append(br)
    for pc,PB in solve(tuple(A),tuple(rhs),len(piv)):
        k=len(PB[0]) if PB else 0
        cc=[(v+sum(x*y for x,y in zip(row,pc)))%1 for v,row in zip(c,B)]
        BB=[tuple(sum(row[h]*PB[h][j] for h in range(len(piv))) for j in range(k)) for row in B]
        yield canon(cc,BB)

def main(name):
    started=time.time();s=family(name);d=s['parameters'];pairs=set()
    e1=[edges(p,s['constants'][0],s['basis'][0]) for p in PAT]
    e2=[edges(p,s['constants'][1],s['basis'][1]) for p in PAT]
    for p in e1:
        for q in e2:pairs.update(joint(p+q,d))
    print(name,'independent pair components',len(pairs),flush=True)
    out=set()
    for num,comp in enumerate(pairs,1):
        for p in PAT:out.update(restrict(comp,p,s))
        if num%1000==0:print(name,'independent progress',num,'/',len(pairs),flush=True)
    maximal=[]
    for dim in sorted({len(c[2]) for c in out},reverse=True):
        layer=[c for c in out if len(c[2])==dim]
        maximal.extend(c for c in layer if not any(contains(big,c) for big in maximal))
    source=json.loads((ROOT/'results'/(name+'.json')).read_text())
    expected={deserialize(r) for r in source['maximal_components']}
    assert set(maximal)==expected,(name,len(maximal),len(expected))
    assert len(out)==source['distinct_final_components']
    r=dict(name=name,pair_components=len(pairs),distinct_final_components=len(out),maximal_components=len(maximal),dimensions=dict(Counter(len(c[2]) for c in maximal)),match=True,seconds=time.time()-started)
    (ROOT/'results'/(name+'_independent.json')).write_text(json.dumps(r,indent=2)+'\n')
    print(name,'INDEPENDENT MATCH',r,flush=True)
if __name__=='__main__':
    for name in sys.argv[1:] or ['F4','F5','F6','F2','F3','F1']:main(name)
