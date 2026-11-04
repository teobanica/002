"""Regression without support symmetry: all 105 oriented C supports.

Uses the unimodular solver, but constructs supports/phase equations afresh.
Independent of the nine-case reduction, not a third arithmetic solver.
"""
from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
import json
from torus_solver import solve_torus

ROOT=Path(__file__).resolve().parent


def key(H):
    zero=(F(0),)*7
    H=(zero,)+tuple(H)
    result=[]
    for b,r,s in permutations(range(3)):
        columns=[((H[r][j]-H[b][j])%1,(H[s][j]-H[b][j])%1) for j in range(7)]
        for x,y in columns:
            result.append(tuple(sorted(((a-x)%1,(c-y)%1) for a,c in columns)))
    return min(result)


def main():
    verified=json.loads((ROOT/'verification.json').read_text())
    result={}
    for kind,p in (('ACC',7),('BCC',5)):
        m=verified['types'][kind]['matrix_modulus']
        targets={key([[F(x,m) for x in r] for r in item['exponents'][1:]]) for item in verified['types'][kind]['representatives']}
        patterns=[]
        for T in combinations(range(7),p):
            rest=tuple(j for j in range(7) if j not in T)
            for powers in permutations(range(1,p)):
                edges=[(T[0],j,F(e,p)) for j,e in zip(T[1:],powers)]
                if rest:edges.append((rest[0],rest[1],F(1,2)))
                patterns.append(edges)
        count=consistent=solutions=0
        point_set=set()
        supports_seen=set()
        for T in combinations(range(7),3):
            rest=set(range(7))-set(T)
            for P in combinations(sorted(rest),2):
                Q=tuple(sorted(rest-set(P)))
                if P>Q:continue
                supports_seen.add((T,P,Q))
                for sign in (1,-1):
                    uc=[F(0),F(1,3),F(2,3),F(0),F(1,2),F(0),F(1,2)]
                    UB=[[0]*4 for _ in range(7)]
                    for j in (3,4):UB[j][0]=1
                    for j in (5,6):UB[j][1]=1
                    vc=[F(0)]*7;VB=[[0]*4 for _ in range(7)]
                    for e,j in enumerate(T):vc[j]=F(sign*e,3)%1
                    for v,pair in ((2,P),(3,Q)):
                        for j in pair:VB[j][v]=1
                        vc[pair[1]]=F(1,2)
                    c=[a-b for a,b in zip(uc,vc)]
                    B=[[a-b for a,b in zip(r,s)] for r,s in zip(UB,VB)]
                    for edges in patterns:
                        A=[[B[j][k]-B[i][k] for k in range(4)] for i,j,_ in edges]
                        rhs=[(t-c[j]+c[i])%1 for i,j,t in edges]
                        points=solve_torus(A,rhs,4)
                        count+=1;consistent+=bool(points);solutions+=len(points)
                        for x,K in points:
                            assert len(K[0])==0
                            H=[]
                            for rc,R in ((uc,UB),(vc,VB)):
                                values=[(rc[j]+sum(R[j][k]*x[k] for k in range(4)))%1 for j in range(7)]
                                H.append(tuple((a-values[0])%1 for a in values))
                            point_set.add(tuple(H))
        assert len(supports_seen)==105
        found={key(H) for H in point_set}
        assert found==targets
        assert (count,consistent,solutions)==((151200,576,576) if kind=='ACC' else (105840,384,768))
        result[kind]={'systems':count,'consistent_systems':consistent,'point_occurrences':solutions,
                      'distinct_normalized_points':len(point_set),'classes':len(found)}
        print('All 105 supports:',kind,result[kind],flush=True)
    (ROOT/'full_support_check.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
