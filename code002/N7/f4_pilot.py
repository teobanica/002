"""Exact cycle-equation enumeration of fourth-row extensions.
Standard library only. Angles are in R/Z; x_0=0 fixes the new row phase.
"""
from itertools import combinations, permutations, product
from fractions import Fraction as Q
from collections import Counter
from pathlib import Path
import json, time

MOD=840 # accommodates roots of orders 24, 7, 5, 3, 2

def patterns(kind):
    out=[]
    def add(groups):
        edges=[]
        for indices, residues, order in groups:
            a=indices[0]
            for j,k in enumerate(indices[1:],1):
                edges.append((k,a,(residues[j]-residues[0])*(MOD//order)))
        out.append(tuple(edges))
    if kind=='A':
        for p in permutations(range(1,7)):
            add([(tuple(range(7)),(0,)+p,7)])
    if kind=='B':
        for pair in combinations(range(7),2):
            five=tuple(j for j in range(7) if j not in pair)
            for p in permutations(range(1,5)):
                add([(five,(0,)+p,5),(pair,(0,1),2)])
    if kind=='C':
        for tri in combinations(range(7),3):
            rem=[j for j in range(7) if j not in tri]
            for k in rem[1:]:
                pair=(rem[0],k); other=tuple(j for j in rem if j not in pair)
                for ori in [(0,1,2),(0,2,1)]:
                    add([(tri,ori,3),(pair,(0,1),2),(other,(0,1),2)])
    return out
PAT={t:patterns(t) for t in 'ABC'}

class State:
    def __init__(self):
        self.p=list(range(7));self.c=[0]*7;self.n=[0]*7
        self.q=None # None = unrestricted parameter; else a tuple of exact angles
    def clone(self):
        s=State();s.p=self.p.copy();s.c=self.c.copy();s.n=self.n.copy();s.q=self.q
        return s
    def find(self,i):
        c=n=0
        while self.p[i]!=i:
            c+=self.c[i];n+=self.n[i];i=self.p[i]
        return i,c%MOD,n
    def close(self,c,n):
        c%=MOD
        if self.q is not None:
            self.q=tuple(t for t in self.q if (Q(c,MOD)+n*t)%1==0)
            return bool(self.q)
        if n==0:return c==0
        if n<0:c,n=-c,-n
        self.q=tuple(sorted({((Q(-c,MOD)+k)/n)%1 for k in range(n)}))
        return True
    def add(self,i,j,c,n):
        a,ca,na=self.find(i);b,cb,nb=self.find(j)
        c=(c-ca+cb)%MOD;n=n-na+nb
        if a==b:return self.close(c,n)
        if a<b:a,b,c,n=b,a,-c,-n
        self.p[a]=b;self.c[a]=c%MOD;self.n[a]=n
        return True
    def apply(self,edges,row):
        for i,j,c in edges:
            if not self.add(i,j,c+row[i][0]-row[j][0],row[i][1]-row[j][1]):return False
        return True
    def keys(self):
        labels=[self.find(i) for i in range(7)]
        if self.q is None:
            yield ('generic',tuple(labels)),self
        else:
            for q in self.q:
                key=('fixed',q,tuple((r,(Q(c,MOD)+n*q)%1) for r,c,n in labels))
                s=self.clone();s.q=(q,)
                yield key,s
    def serialize(self):
        labels=[self.find(i) for i in range(7)]
        roots=sorted({v[0] for v in labels});free=roots[1:]
        q=None if self.q is None else self.q[0]
        return dict(parameter=None if q is None else str(q),
             free_components=free, row=[dict(component=r,constant=str(Q(c,MOD) if q is None else (Q(c,MOD)+n*q)%1),slope=n if q is None else 0) for r,c,n in labels])

def seed(name):
    if name=='E9':
        return [[(0,0)]*7,[(x*35,0) for x in [0,8,16,2,14,6,18]],[(x*35,0) for x in [0,4,2,10,16,18,12]]]
    if name=='F4':
        return [[(0,0)]*7,[(x*140,int(j>=3)) for j,x in enumerate([0,2,4,2,5,0,3])],[(x*140,int(j>=3)) for j,x in enumerate([0,2,4,0,2,3,5])]]
    raise ValueError(name)

# For a CCC seed, allowed triangle types imply at most one of the three
# new edges can be non-C. Types ACC, BCC, CCC are all explicitly enumerated.
TYPES=['CCC','ACC','CAC','CCA','BCC','CBC','CCB']

def enumerate_seed(name):
    rows=seed(name); final={};stats=[];started=time.time()
    for types in TYPES:
        # Put the two C edges first, to keep the intermediate list small.
        order=sorted(range(3),key=lambda j:types[j]!='C')
        unique={};counts=[0,0,0]
        for i,p in enumerate(PAT[types[order[0]]]):
            s=State();assert s.apply(p,rows[order[0]])
            for j,t in enumerate(PAT[types[order[1]]]):
                counts[0]+=1
                z=s.clone()
                if not z.apply(t,rows[order[1]]):continue
                counts[1]+=1
                for key,z1 in z.keys():unique.setdefault(key,(z1,[(order[0],i),(order[1],j)]))
        before=len(final);success=0
        for s,wit in unique.values():
            for k,t in enumerate(PAT[types[order[2]]]):
                counts[2]+=1;z=s.clone()
                if not z.apply(t,rows[order[2]]):continue
                success+=1
                for key,z1 in z.keys():
                    if key not in final:
                        d=z1.serialize();d['types']=types;d['witness']=sorted(wit+[(order[2],k)])
                        final[key]=d
        stat=dict(types=types,pair_systems=counts[0],consistent_pair_systems=counts[1],distinct_pair_loci=len(unique),third_systems=counts[2],consistent_third_systems=success,new_loci=len(final)-before)
        stats.append(stat);print(name,stat,flush=True)
    result=dict(seed=name,statistics=stats,loci=list(final.values()),seconds=time.time()-started)
    (Path(__file__).parent/(name+'_results.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(name,'TOTAL',len(final),'dimensions',Counter((d['parameter'] is None,len(d['free_components'])) for d in final.values()),flush=True)
    return result

if __name__=='__main__':
    import sys
    for name in sys.argv[1:] or ['E9','F4']:enumerate_seed(name)
