"""Independent exact witness, locus cross-check and finite-graph verification."""
import json
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from math import lcm
from torus import *
R=Path(__file__).resolve().parents[1]
def profile(v,primes):
    c=Counter(x%1 for x in v)
    def rec(c,pp):
        if not pp:return not c
        x=min(c)
        for p in set(pp):
            cc=c.copy()
            for k in range(p):
                y=(x+F(k,p))%1
                if not cc[y]:break
                cc[y]-=1
                if not cc[y]:del cc[y]
            else:
                pr=list(pp);pr.remove(p)
                if rec(cc,pr):return True
        return False
    return rec(c,list(primes))
def regular(r,s):
    v=[x-y for x,y in zip(r,s)]
    return any(profile(v,p) for p in [(2,2,2,2),(2,3,3),(3,5)])
def maxclique(adj):
    best=[]
    def f(chosen,pool):
        nonlocal best
        if len(chosen)>len(best):best=chosen
        while len(chosen)+len(pool)>len(best):
            v=pool.pop();f(chosen+[v],[u for u in pool if u in adj[v]])
    f([],list(range(len(adj))));return best
def symbolic_A(r,s,mod):
    c=Counter(((a[0]-b[0])%mod,)+tuple(x-y for x,y in zip(a[1:],b[1:])) for a,b in zip(r,s))
    return all(n==c[((e[0]+mod//2)%mod,)+e[1:]] for e,n in c.items())
L=json.loads((R/'outputs/class21_line_parameters.json').read_text());S=json.loads((R/'outputs/class21_line_scan.json').read_text())
lm=int((R/'outputs/line_modulus.txt').read_text())
for l,s in zip(L,S):
    if s['k5']:
        H=[[[int(F(*z[0])*lm),z[1]] for z in r] for r in l['rows']]+[s['rows'][j] for j in s['k5']]
        assert all(symbolic_A(r,s,lm) for r,s in combinations(H,2))
print('PASS: 16 symbolic square witnesses (8 circles in each component).')
EV=json.loads((R/'outputs/class21_A_events.json').read_text());residue=json.loads((R/'outputs/final_residue.json').read_text())
C=[c for c in json.loads((R/'outputs/components.json').read_text()) if c['class_id']==21]
summary=[];mod=1
for br in range(2):
    D=json.loads((R/f'outputs/class21_full{br}.json').read_text())
    a=D['profiles']['AAA'];assert set(map(str,a['lines']))==set(map(str,EV[br]['lines']));assert set(map(str,a['points']))==set(map(str,EV[br]['points']))
    assert not any('C' in t for t in D['profiles'])
    assert all(not r['global'] and not r['lines'] for t,r in D['profiles'].items() if t!='AAA')
    cc=Counter();rowcounts=Counter();fib={decode_point(f['point']):f for f in D['fibers']}
    for i,ep in enumerate(residue[br]['points']):
        p=decode_point(ep);f=fib[p];assert f['disconnected']==0
        rows=[tuple(F(*q) for q in r) for r in f['rows']]
        for r in rows:
            for q in r:mod=lcm(mod,q.denominator)
        H=[[(F(z,30)+a*p[0]+b*p[1])%1 for z,a,b in r] for r in C[br]['rows']]
        assert all(regular(r,s) for r in rows for s in H)
        adj=[set() for _ in rows]
        for j,k in combinations(range(len(rows)),2):
            if regular(rows[j],rows[k]):adj[j].add(k);adj[k].add(j)
        cl=maxclique(adj);assert len(cl)<5,(br,p,cl)
        cc[len(cl)]+=1;rowcounts[len(rows)]+=1
    summary.append(dict(branch=br,points=len(fib),row_count_distribution=dict(rowcounts),clique_number_distribution=dict(cc)))
    print('PASS: branch',br,'points',len(fib),'candidate counts',dict(rowcounts),'clique numbers',dict(cc))
(R/'outputs/final_verification.json').write_text(json.dumps(dict(components=summary,row_root_order_lcm=mod),indent=2))
print('All finite candidates lie in mu_'+str(mod)+'.')
for br in range(2):
    D=json.loads((R/f'outputs/class21_full{br}.json').read_text());finite=json.loads((R/f'outputs/finite{br}.json').read_text())
    fibers={decode_point(f['point']):{tuple(F(*q) for q in r) for r in f['rows']} for f in D['fibers']}
    assert len(finite)==64
    for rr in finite:
        p=decode_point(residue[br]['points'][rr['id']]);actual={tuple(F(q,rr['modulus']) for q in r) for r in rr['rows']}
        assert actual==fibers[p],(br,p,len(actual),len(fibers[p]))
print('PASS: independent finite-root enumeration reproduces all 128 candidate lists exactly.')
import csv
p=R/'outputs/AAA_status.csv';table=list(csv.DictReader(p.open()));keys=list(table[0])
for row in table:
    if row['class']=='21':row['general_completion_status']='exact eight-circle iff locus proved in this run; see class-21 theorem'
with p.open('w') as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(table)
