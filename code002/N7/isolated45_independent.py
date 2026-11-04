"""Independent verification WITHOUT any three-row type exclusions.
Handles every new A edge directly, then exhausts all B/C choices on all
three new edges. Generates cycle patterns independently from 7! permutations
and uses graph traversal, not weighted merging. All phases remain unrestricted.
"""
from pathlib import Path
from itertools import permutations,combinations
from collections import Counter,deque
import json,time
from seeds import MOD,seeds
BASE=Path(__file__).resolve().parent

def make_patterns(t):
    out=set();blocks={'B':[(0,5),(5,7)],'C':[(0,3),(3,5),(5,7)]}[t]
    for p in permutations(range(7)):
        edges=[]
        for lo,hi in blocks:
            b=p[lo:hi];m=hi-lo;root=min(b);anchor=b.index(root)
            edges.extend((j,root,((k-anchor)%m)*(MOD//m)) for k,j in enumerate(b) if j!=root)
        out.add(tuple(sorted(edges)))
    return sorted(out)
PATS={t:make_patterns(t) for t in 'BC'}
assert len(PATS['B'])==504 and len(PATS['C'])==210
PATS=PATS['B']+PATS['C']

def regular(v):
    count=Counter(x%MOD for x in v)
    if len(count)==7 and all(count[(v[0]+k*(MOD//7))%MOD]==1 for k in range(7)):return 'A'
    for m in [5,3]:
        for root in count:
            cycle=[(root+k*(MOD//m))%MOD for k in range(m)]
            if not all(count[x] for x in cycle):continue
            rest=count.copy();rest.subtract(cycle)
            if all(rest[x]==rest[(x+MOD//2)%MOD] for x in rest):return 'B' if m==5 else 'C'
    return None

def edges(pattern,row):return [(i,j,(c+row[i]-row[j])%MOD) for i,j,c in pattern]

def solve(eq):
    adj=[[] for _ in range(7)]
    for i,j,c in eq:adj[j].append((i,c));adj[i].append((j,-c))
    labels=[None]*7
    for root in range(7):
        if labels[root] is not None:continue
        labels[root]=(root,0);todo=[root]
        while todo:
            j=todo.pop();c0=labels[j][1]
            for i,c in adj[j]:
                c=(c+c0)%MOD
                if labels[i] is None:labels[i]=(root,c);todo.append(i)
                elif labels[i][1]!=c:return None
    return tuple(labels)

def main():
    report=[]
    tails=list(permutations(range(1,7)))
    for seed in seeds():
        started=time.time();rows=seed['rows'];found=set();a_candidates=set()
        for old in rows:
            for p in tails:
                row=tuple((old[j]-old[0]+r*(MOD//7))%MOD for j,r in enumerate((0,)+p))
                a_candidates.add(row)
        a_found=0
        for row in a_candidates:
            if all(regular([x-y for x,y in zip(row,old)]) for old in rows):
                found.add(tuple((0,x) for x in row));a_found+=1
        eds=[[edges(p,row) for p in PATS] for row in rows]
        pairloci={};consistentpairs=0
        for p in eds[0]:
            for q in eds[1]:
                key=solve(p+q)
                if key is not None:
                    consistentpairs+=1
                    pairloci.setdefault(key,[(i,r,c) for i,(r,c) in enumerate(key) if i!=r])
        consistentthird=0
        for p in pairloci.values():
            for q in eds[2]:
                key=solve(p+q)
                if key is not None:consistentthird+=1;found.add(key)
        primary=json.loads((BASE/'results'/(seed['label']+'.json')).read_text())
        expected={tuple(zip(d['components'],d['constants'])) for d in primary['loci']}
        assert found==expected,(seed['label'],found-expected,expected-found)
        r=dict(label=seed['label'],a_candidates=len(a_candidates),a_survivors=a_found,
               bc_pair_systems=len(PATS)**2,consistent_bc_pairs=consistentpairs,
               bc_pair_loci=len(pairloci),bc_third_systems=len(pairloci)*len(PATS),
               consistent_bc_thirds=consistentthird,solution_loci=len(found),
               positive_dimensional_loci=sum(len({r for r,c in k})>1 for k in found),
               matches_primary=True,seconds=round(time.time()-started,3))
        report.append(r)
        (BASE/'independent_results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(seed['label'],'matched',len(found),'loci;',r['bc_pair_loci'],'pair loci;',r['seconds'],'sec',flush=True)
    print('ALL 45 AGREE, WITHOUT USING THREE-ROW TYPE EXCLUSIONS.',flush=True)

if __name__=='__main__':main()
