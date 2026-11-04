"""Independent direct check of CCC validity and ordinary monomial equivalence.

No cycle-orbit pruning routine is imported. All 3!*7! transformations
are visited for every representative; only explicit necessary equalities
are used to reject a transformation early.
"""
from itertools import permutations,combinations
from pathlib import Path
from collections import Counter
import json,time
from families import templates
ROOT=Path(__file__).parent

def blocks(c,B):
    for T in combinations(range(7),3):
        if not B[T[0]]==B[T[1]]==B[T[2]]:continue
        for s in (1,-1):
            if any((c[T[k]]-c[T[0]]-8*s*k)%24 for k in (1,2)):continue
            a,b,d,e=[j for j in range(7) if j not in T]
            for P,Q in [((a,b),(d,e)),((a,d),(b,e)),((a,e),(b,d))]:
                if all(B[i]==B[j] and (c[j]-c[i])%24==12 for i,j in (P,Q)):
                    return {'triple':T,'orientation':s,'pairs':[P,Q]}
    return None

def fit(target,C,D):
    K,V,P=target
    for i in range(14):
        if (C[i]-K[i]-sum(V[i][h]*C[p] for h,p in enumerate(P)))%24:return False
        for j in range(len(D[0])):
            if D[i][j]!=sum(V[i][h]*D[p][j] for h,p in enumerate(P)):return False
    return True

def run():
    data=json.loads((ROOT/'coverage.json').read_text())
    targets=[]
    for x in data['normalized_targets']:
        f=x['family'];targets.append((x['name'],(tuple(f['constant_24ths']),tuple(tuple(r) for r in f['basis']),tuple(f['pivot']))))
    report={'checks':[],'validity':[]}
    for name,f in templates():
        c,B,P=f;d=len(P)
        quotient=(tuple((c[j]-c[j+7])%24 for j in range(7)),
                  tuple(tuple(x-y for x,y in zip(B[j],B[j+7])) for j in range(7)))
        cert=[blocks(c[:7],B[:7]),blocks(c[7:],B[7:]),blocks(*quotient)]
        assert all(cert)
        report['validity'].append({'name':name,'cycles':cert})
        R=[(0,)*7,c[:7],c[7:]];S=[((0,)*d,)*7,B[:7],B[7:]]
        visited=0;survived=0
        for base,a,b in permutations(range(3)):
            U=[(R[a][j]-R[base][j])%24 for j in range(7)]
            V=[(R[b][j]-R[base][j])%24 for j in range(7)]
            E=[tuple(S[a][j][k]-S[base][j][k] for k in range(d)) for j in range(7)]
            F=[tuple(S[b][j][k]-S[base][j][k] for k in range(d)) for j in range(7)]
            for anchor in range(7):
                UC=[(x-U[anchor])%24 for x in U];VC=[(x-V[anchor])%24 for x in V]
                UB=[tuple(x-y for x,y in zip(r,E[anchor])) for r in E]
                VB=[tuple(x-y for x,y in zip(r,F[anchor])) for r in F]
                for tail in permutations([j for j in range(7) if j!=anchor]):
                    visited+=1
                    if UC[tail[0]]!=8 or any(UB[tail[0]]):continue
                    if UC[tail[1]]!=16 or any(UB[tail[1]]):continue
                    if any((UC[tail[j+1]]-UC[tail[j]])%24!=12 or UB[tail[j+1]]!=UB[tail[j]]
                           for j in (2,4)):continue
                    survived+=1
                    cols=(anchor,)+tail
                    C=tuple(UC[j] for j in cols)+tuple(VC[j] for j in cols)
                    D=tuple(UB[j] for j in cols)+tuple(VB[j] for j in cols)
                    for other,g in targets:
                        if other!=name and len(g[2])>=d:
                            assert not fit(g,C,D),(name,'contained in',other)
        assert visited==30240
        report['checks'].append({'name':name,'dimension':d,'permutations_visited':visited,
                                 'normalizations_surviving':survived})
        print(name,'valid CCC; all',visited,'transformations checked;',survived,'survive necessary conditions.',flush=True)
    report['total_permutations']=sum(x['permutations_visited'] for x in report['checks'])
    (ROOT/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('TOTAL',report['total_permutations'],'ordinary monomial transformations checked.',flush=True)

if __name__=='__main__':
    t=time.perf_counter();run();print('Elapsed',round(time.perf_counter()-t,2),'seconds.')

