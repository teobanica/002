"""Exhaustive fourth-row loci over each affine seed, exact rational phases."""
from fractions import Fraction as Q
from itertools import combinations
from functools import lru_cache
from collections import Counter
from pathlib import Path
import json,time,sys
from torus_solver import solve_torus
from families import family
ROOT=Path(__file__).resolve().parent

def patterns():
    out=[]
    for tri in combinations(range(7),3):
        rest=[i for i in range(7) if i not in tri]
        for j in rest[1:]:
            pair=(rest[0],j);other=tuple(i for i in rest if i not in pair)
            for ori in [1,-1]:
                out.append(((tri[1],tri[0],(2*ori)%6),(tri[2],tri[0],(4*ori)%6),(pair[1],pair[0],3),(other[1],other[0],3)))
    return out
PAT=patterns()

def det(M):
    if not M:return 1
    if len(M)==1:return M[0][0]
    return sum((-1)**j*M[0][j]*det([r[:j]+r[j+1:] for r in M[1:]]) for j in range(len(M)))

@lru_cache(None)
def chart(B):
    d=len(B[0])
    if d==0:return (),tuple(() for _ in B)
    piv=next((p for p in combinations(range(len(B)),d) if abs(det([list(B[i]) for i in p]))==1),None)
    if piv is None:raise RuntimeError(('No unimodular coordinate chart',B))
    M=[list(B[i]) for i in piv];D=det(M)
    inv=[[(-1)**(i+j)*det([r[:i]+r[i+1:] for k,r in enumerate(M) if k!=j])//D for j in range(d)] for i in range(d)]
    C=tuple(tuple(sum(row[k]*inv[k][j] for k in range(d)) for j in range(d)) for row in B)
    return piv,C

def canonical(c,B):
    B=tuple(tuple(v) for v in B);piv,C=chart(B)
    cc=tuple((v-sum(row[k]*c[piv[k]] for k in range(len(piv))))%1 for v,row in zip(c,C))
    return cc,C,piv

def contains(big,small):
    bc,B,p=big;sc,S,_=small
    if len(S[0])>len(p):return False
    for i,row in enumerate(B):
        if (sc[i]-bc[i]-sum(v*sc[j] for v,j in zip(row,p)))%1:return False
        if any(v!=sum(row[k]*S[j][h] for k,j in enumerate(p)) for h,v in enumerate(S[i])):return False
    return True

@lru_cache(None)
def solve_cached(A,b,n):
    if n==0:return (((),()),) if all(x%1==0 for x in b) else ()
    return tuple(solve_torus(A,b,n))

class Graph:
    def __init__(self,d):
        self.d=d;self.p=list(range(7));self.c=[0]*7;self.B=[(0,)*d for _ in range(7)];self.eq={}
    def clone(self):
        g=Graph(self.d);g.p=self.p.copy();g.c=self.c.copy();g.B=self.B.copy();g.eq=self.eq.copy();return g
    def find(self,i):
        c=0;b=[0]*self.d
        while self.p[i]!=i:
            c+=self.c[i];b=[x+y for x,y in zip(b,self.B[i])];i=self.p[i]
        return i,c%6,tuple(b)
    def add(self,i,j,c,b):
        ri,ci,bi=self.find(i);rj,cj,bj=self.find(j)
        c=(c-ci+cj)%6;b=tuple(x-y+z for x,y,z in zip(b,bi,bj))
        if ri==rj:
            if not any(b):return c==0
            if next(v for v in b if v)<0:c=(-c)%6;b=tuple(-v for v in b)
            if b in self.eq:return self.eq[b]==c
            self.eq[b]=c;return True
        if ri<rj:ri,rj,c,b=rj,ri,-c,tuple(-v for v in b)
        self.p[ri]=rj;self.c[ri]=c%6;self.B[ri]=b;return True
    def apply(self,pat,rc,rb):
        for i,j,c in pat:
            if not self.add(i,j,c+rc[i]-rc[j],tuple(x-y for x,y in zip(rb[i],rb[j]))):return False
        return True
    def joint_components(self):
        eq=sorted(self.eq.items());A=tuple(b for b,c in eq);rhs=tuple(Q(-c,6)%1 for b,c in eq)
        labs=[self.find(j) for j in range(7)];free=sorted({r for r,c,b in labs}-{0})
        for pc,PB in solve_cached(A,rhs,self.d):
            k=len(PB[0]);C=list(pc);B=[row+(0,)*len(free) for row in PB]
            for r,c,b in labs:
                C.append((Q(c,6)+sum(x*y for x,y in zip(b,pc)))%1)
                B.append(tuple(sum(b[i]*PB[i][j] for i in range(self.d)) for j in range(k))+tuple(int(r==s) for s in free))
            yield canonical(C,B)

def restrict(comp,pat,seed):
    cc,B,piv=comp;d=seed['parameters'];A=[];rhs=[]
    for i,j,c in pat:
        coeff=tuple(seed['basis'][2][i][h]-seed['basis'][2][j][h] for h in range(d))
        ac=cc[d+i]-cc[d+j]-sum(coeff[h]*cc[h] for h in range(d))
        ar=tuple(B[d+i][k]-B[d+j][k]-sum(coeff[h]*B[h][k] for h in range(d)) for k in range(len(piv)))
        target=(Q(c+seed['constants'][2][i]-seed['constants'][2][j],6)-ac)%1
        if not any(ar):
            if target:return
        else:A.append(ar);rhs.append(target)
    for c,D in solve_cached(tuple(A),tuple(rhs),len(piv)):
        k=len(D[0]) if D else 0
        outc=[(v+sum(x*y for x,y in zip(row,c)))%1 for v,row in zip(cc,B)]
        outB=[tuple(sum(row[h]*D[h][j] for h in range(len(piv))) for j in range(k)) for row in B]
        yield canonical(outc,outB)

def serialize(comp):
    c,B,p=comp;return dict(constant=[str(v) for v in c],basis=B,pivot=p,dimension=len(p))
def deserialize(d):return tuple(Q(v) for v in d['constant']),tuple(tuple(v) for v in d['basis']),tuple(d['pivot'])

def main(name):
    start=time.time();seed=family(name);d=seed['parameters'];pairs={};consistent=0
    for i,p in enumerate(PAT):
        g=Graph(d);assert g.apply(p,seed['constants'][0],seed['basis'][0])
        for j,q in enumerate(PAT):
            z=g.clone()
            if not z.apply(q,seed['constants'][1],seed['basis'][1]):continue
            cs=list(z.joint_components())
            if cs:consistent+=1
            for comp in cs:pairs.setdefault(comp,(i,j))
    print(name,'pair systems',210**2,'consistent',consistent,'distinct components',len(pairs),flush=True)
    final={};surviving_systems=0
    for num,(comp,ij) in enumerate(pairs.items(),1):
        for k,pat in enumerate(PAT):
            sol=list(restrict(comp,pat,seed))
            if sol:surviving_systems+=1
            for out in sol:final.setdefault(out,ij+(k,))
        if num%500==0:print(name,'processed pairs',num,'/',len(pairs),'final components',len(final),flush=True)
    maximal=[]
    for dim in sorted({len(c[2]) for c in final},reverse=True):
        layer=[c for c in final if len(c[2])==dim]
        kept=[c for c in layer if not any(contains(big,c) for big in maximal)]
        maximal+=kept
    maximal.sort(key=lambda v:(-len(v[2]),v))
    records=[dict(serialize(c),cycle_patterns=final[c]) for c in maximal]
    result=dict(family=seed,pair_systems=210**2,consistent_pair_systems=consistent,pair_components=len(pairs),
                restricted_systems=len(pairs)*210,consistent_restricted_systems=surviving_systems,
                distinct_final_components=len(final),raw_dimensions=dict(Counter(len(c[2]) for c in final)),
                maximal_components=records,maximal_dimensions=dict(Counter(len(c[2]) for c in maximal)),seconds=time.time()-start)
    (ROOT/'results'/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(name,'DONE',len(maximal),'maximal',result['maximal_dimensions'],'raw',len(final),'seconds',result['seconds'],flush=True)
    return result
if __name__=='__main__':
    for name in sys.argv[1:] or ['F4','F5','F6','F2','F3','F1']:main(name)
