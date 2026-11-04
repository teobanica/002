"""Independent graph-traversal solver; does not use pilot.State.
Regenerates cycle patterns from column permutations. Repeats the exhaustive
extension enumeration and compares the resulting loci with the main solver.
"""
from itertools import permutations
from fractions import Fraction as Q
from collections import deque
from pathlib import Path
import json,time
from pilot import MOD, seed, TYPES

def make_patterns(t):
    out=set()
    blocks={'A':[(0,7)],'B':[(0,5),(5,7)],'C':[(0,3),(3,5),(5,7)]}[t]
    for p in permutations(range(7)):
        groups=[]
        for lo,hi in blocks:
            v=p[lo:hi];m=hi-lo;r=min(v);offset=v.index(r)
            entries=sorted((j,((k-offset)%m)*(MOD//m)) for k,j in enumerate(v))
            groups.append(tuple((j,r,c) for j,c in entries if j!=r))
        out.add(tuple(sorted(e for g in groups for e in g)))
    return sorted(out)
PAT={t:make_patterns(t) for t in 'ABC'}
assert {t:len(p) for t,p in PAT.items()}=={'A':720,'B':504,'C':210}

def edge_list(pattern,row):
    return [(i,j,(c+row[i][0]-row[j][0])%MOD,row[i][1]-row[j][1]) for i,j,c in pattern]

def solve(edges,qfixed=None):
    adj=[[] for _ in range(7)]
    for i,j,c,n in edges:
        adj[j].append((i,c,n));adj[i].append((j,-c,-n))
    labels=[None]*7;relations=set()
    for root in range(7):
        if labels[root] is not None:continue
        labels[root]=(root,0,0);todo=deque([root])
        while todo:
            j=todo.popleft();_,c0,n0=labels[j]
            for i,c,n in adj[j]:
                c=(c+c0)%MOD;n+=n0
                if labels[i] is None:
                    labels[i]=(root,c,n);todo.append(i)
                else:
                    _,c1,n1=labels[i];dc=(c-c1)%MOD;dn=n-n1
                    if dn==0:
                        if dc:return []
                    else:relations.add((dc,dn))
    qs=None if qfixed is None else {qfixed}
    for c,n in relations:
        if n<0:c,n=-c,-n
        possible={((Q(-c,MOD)+k)/n)%1 for k in range(n)}
        qs=possible if qs is None else qs&possible
        if not qs:return []
    compressed=[(i,r,c,n) for i,(r,c,n) in enumerate(labels) if i!=r]
    if qs is None:return [(('generic',tuple(labels)),compressed,None)]
    return [(('fixed',q,tuple((r,(Q(c,MOD)+n*q)%1) for r,c,n in labels)),compressed,q) for q in sorted(qs)]

def main(name):
    rows=seed(name);allkeys=set();stats=[]
    for typ in TYPES:
        order=sorted(range(3),key=lambda i:typ[i]!='C')
        eds=[[edge_list(p,rows[i]) for p in PAT[typ[i]]] for i in order]
        pairs={}
        for e in eds[0]:
            for f in eds[1]:
                for k,comp,q in solve(e+f):pairs[k]=(comp,q)
        survivors=0
        for e,q in pairs.values():
            for f in eds[2]:
                solutions=solve(e+f,q)
                survivors+=bool(solutions)
                allkeys.update(k for k,_,_ in solutions)
        stats.append(dict(types=typ,pair_loci=len(pairs),surviving_third_systems=survivors))
        print(name,stats[-1],flush=True)
    expected=set()
    data=json.loads((Path(__file__).parent/(name+'_results.json')).read_text())
    for d in data['loci']:
        if d['parameter'] is None:
            expected.add(('generic',tuple((v['component'],int(Q(v['constant'])*MOD),v['slope']) for v in d['row'])))
        else:
            expected.add(('fixed',Q(d['parameter']),tuple((v['component'],Q(v['constant'])) for v in d['row'])))
    assert allkeys==expected,(name,len(allkeys),len(expected),allkeys-expected,expected-allkeys)
    print(name,'independent match:',len(allkeys),'loci',flush=True)
    return dict(seed=name,statistics=stats,loci=len(allkeys),match=True)

if __name__=='__main__':
    out=[main(n) for n in ['E9','F4']]
    (Path(__file__).parent/'independent_results.json').write_text(json.dumps(out,indent=2)+'\n')
