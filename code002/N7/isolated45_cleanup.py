"""Complete regular fourth-row extension and compatibility graphs for 45 seeds.
Exact modular graph equations; disconnected graphs retain free circle phases.
"""
from itertools import combinations,permutations,product
from collections import Counter
from pathlib import Path
import json,time,sys
from seeds import MOD,seeds
BASE=Path(__file__).resolve().parent
ALLOWED={'AAA','ACC','BBB','BCC','CCC'}

def patterns(t):
    out=[]
    def add(groups):
        edges=[]
        for b,res,m in groups:
            for j,r in zip(b[1:],res[1:]):edges.append((j,b[0],r*(MOD//m)))
        out.append(tuple(edges))
    if t=='A':
        for p in permutations(range(1,7)):add([(tuple(range(7)),(0,)+p,7)])
    if t=='B':
        for pair in combinations(range(7),2):
            rest=tuple(j for j in range(7) if j not in pair)
            for p in permutations(range(1,5)):add([(rest,(0,)+p,5),(pair,(0,1),2)])
    if t=='C':
        for tri in combinations(range(7),3):
            rest=[j for j in range(7) if j not in tri]
            for j in rest[1:]:
                pair=(rest[0],j);other=tuple(k for k in rest if k not in pair)
                for p in [(0,1,2),(0,2,1)]:add([(tri,p,3),(pair,(0,1),2),(other,(0,1),2)])
    return out
PAT={t:patterns(t) for t in 'ABC'}

def block(row,b,m):return {(row[j]-row[b[0]])%MOD for j in b}==set(range(0,MOD,MOD//m))

def direct_type(row):
    if block(row,tuple(range(7)),7):return 'A',[[j+1 for j in range(7)]]
    for pair in combinations(range(7),2):
        if not block(row,pair,2):continue
        rest=tuple(j for j in range(7) if j not in pair)
        if block(row,rest,5):return 'B',[[j+1 for j in rest],[j+1 for j in pair]]
    for tri in combinations(range(7),3):
        if not block(row,tri,3):continue
        rest=[j for j in range(7) if j not in tri]
        for j in rest[1:]:
            pair=(rest[0],j);other=tuple(k for k in rest if k not in pair)
            if block(row,pair,2) and block(row,other,2):return 'C',[[k+1 for k in b] for b in (tri,pair,other)]
    return None,None

def quotient(a,b):return [(x-y)%MOD for x,y in zip(a,b)]

def assignments(rows):
    old={(i,j):direct_type(quotient(rows[i],rows[j]))[0] for i,j in combinations(range(3),2)}
    assert None not in old.values()
    return [''.join(t) for t in product('ABC',repeat=3) if all(''.join(sorted([old[i,j],t[i],t[j]])) in ALLOWED for i,j in old)],old

class Graph:
    def __init__(self):self.p=list(range(7));self.c=[0]*7
    def clone(self):
        q=Graph();q.p=self.p.copy();q.c=self.c.copy();return q
    def find(self,i):
        c=0
        while i!=self.p[i]:c+=self.c[i];i=self.p[i]
        return i,c%MOD
    def add(self,i,j,c):
        a,x=self.find(i);b,y=self.find(j);c=(c-x+y)%MOD
        if a==b:return c==0
        if a<b:a,b,c=b,a,-c
        self.p[a]=b;self.c[a]=c%MOD;return True
    def apply(self,edges,row):
        for i,j,c in edges:
            if not self.add(i,j,c+row[i]-row[j]):return False
        return True
    def key(self):return tuple(self.find(j) for j in range(7))

def enumerate_seed(seed):
    rows=seed['rows'];allowed,old=assignments(rows)
    assert ''.join(sorted(old.values()))==seed['kind']
    sols={};stats=[];cache={}
    for types in allowed:
        # Start with a seven-cycle if present: it fixes the new row entirely.
        if 'A' in types:
            a=types.index('A');tested=0;passed=0
            for pat in PAT['A']:
                s=Graph();assert s.apply(pat,rows[a]);key=s.key();assert {r for r,c in key}=={0}
                row=[c for r,c in key];tested+=1
                if all(direct_type(quotient(row,rows[i]))[0]==types[i] for i in range(3)):
                    passed+=1;sols.setdefault(key,dict(types=types,origin='A anchor'))
            stats.append(dict(types=types,anchor_row=a+1,tested=tested,consistent=passed));continue
        order=sorted(range(3),key=lambda i:len(PAT[types[i]]))
        a,b,c=order;cachekey=(a,types[a],b,types[b])
        if cachekey not in cache:
            unique={};tested=consistent=0
            for p in PAT[types[a]]:
                s=Graph();assert s.apply(p,rows[a])
                for q in PAT[types[b]]:
                    tested+=1;z=s.clone()
                    if z.apply(q,rows[b]):consistent+=1;unique.setdefault(z.key(),z)
            cache[cachekey]=(unique,tested,consistent)
        unique,tested,consistent=cache[cachekey];third=ok=0
        for s in unique.values():
            for p in PAT[types[c]]:
                third+=1;z=s.clone()
                if z.apply(p,rows[c]):
                    ok+=1;sols.setdefault(z.key(),dict(types=types,origin='three cycle systems'))
        stats.append(dict(types=types,pair_rows=[a+1,b+1],pair_systems=tested,consistent_pair_systems=consistent,pair_loci=len(unique),third_systems=third,consistent_third_systems=ok))
    loci=[]
    for key,meta in sorted(sols.items()):
        roots=sorted({r for r,c in key})
        loci.append(dict(components=[r for r,c in key],constants=[c for r,c in key],dimension=len(roots)-1,**meta))
    return dict(seed=seed,allowed_new_types=allowed,statistics=stats,loci=loci)

def compatible_graph(result):
    if any(d['dimension'] for d in result['loci']):return dict(status='continuous_loci_need_separate_analysis')
    rows=[d['constants'] for d in result['loci']];edges=[];cert=[]
    for i,j in combinations(range(len(rows)),2):
        typ,blocks=direct_type(quotient(rows[i],rows[j]))
        if typ:edges.append((i,j));cert.append(dict(pair=[i+1,j+1],type=typ,cycles=blocks))
    es=set(edges);cliques={0:[()]}
    for k in range(1,5):
        cliques[k]=[c+(j,) for c in cliques[k-1] for j in range(c[-1]+1 if c else 0,len(rows)) if all((i,j) in es for i in c)]
    return dict(status='finite',edges=[[i+1,j+1] for i,j in edges],edge_certificates=cert,
                extension_counts={str(k+3):len(v) for k,v in cliques.items() if k},
                square_completions=[[i+1 for i in c] for c in cliques[4]],
                maximum_rows=3+max(k for k,v in cliques.items() if v))

def main(labels=None):
    selected=[s for s in seeds() if not labels or s['label'] in labels]
    for seed in selected:
        start=time.time();r=enumerate_seed(seed);r['compatibility']=compatible_graph(r);r['seconds']=round(time.time()-start,3)
        (BASE/'results'/(seed['label']+'.json')).write_text(json.dumps(r,indent=2)+'\n')
        print(seed['label'],seed['kind'],'rows:',len(r['loci']),'dimensions:',dict(Counter(d['dimension'] for d in r['loci'])),r['compatibility'].get('extension_counts'),r['seconds'],'sec',flush=True)

if __name__=='__main__':main(sys.argv[1:])
