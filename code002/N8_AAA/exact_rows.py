from fractions import Fraction as F
from collections import Counter
from itertools import combinations

def profile(row,parts):
    """Cyclotomic cycle decomposition, retaining formal monomial exponents."""
    c=Counter((F(e[0])%1,)+tuple(e[1:]) for e in row)
    def rec(c,ps):
        if not ps:return not c
        x=min(c)
        for p in set(ps):
            cc=c.copy();ok=True
            for j in range(p):
                y=((x[0]+F(j,p))%1,)+x[1:]
                if not cc[y]:ok=False;break
                cc[y]-=1
                if not cc[y]:del cc[y]
            if ok:
                pp=list(ps);pp.remove(p)
                if rec(cc,pp):return True
        return False
    return rec(c,parts)
def profiles(row):
    return ''.join(t for t,ps in [('A',[2]*4),('B',[2,3,3]),('C',[3,5])] if profile(row,ps))
def quotient(a,b):return [tuple([((x[0]-y[0])%1)]+[u-v for u,v in zip(x[1:],y[1:])]) for x,y in zip(a,b)]
def clique(adj,need):
    def rec(pool,sol):
        if len(sol)==need:return sol
        while len(pool)+len(sol)>=need:
            v=pool[-1];pool=pool[:-1];r=rec([u for u in pool if adj[v][u]],sol+[v])
            if r:return r
        return []
    return rec(list(range(len(adj))),[])
def graph(rows):
    n=len(rows);adj=[[False]*n for _ in range(n)]
    for i,j in combinations(range(n),2):adj[i][j]=adj[j][i]=bool(profiles(quotient(rows[i],rows[j])))
    k=clique(adj,5)
    maxc=5 if k else next((z for z in range(min(4,n),0,-1) if clique(adj,z)),0)
    return k,maxc,sum(map(sum,adj))//2
