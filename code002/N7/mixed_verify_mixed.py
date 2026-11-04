"""Independent entrywise checks of the literal ACC/BCC normal forms.

No imports from either solver or the support enumerator. Exhaustively tests
all 6*7! row/column permutations for every one of the 56 normalized points.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
ACC_REPS = ((1,3,4,6),(1,5,2,6),(1,6,2,5),(1,6,4,3),
            (2,3,4,5),(2,5,1,6),(2,5,4,3),(2,6,1,5),
            (4,3,1,6),(4,3,2,5),(4,5,2,3),(4,6,1,3))
S = ((1,2),(1,3),(2,1),(3,1))


def M(q, modulus):
    x,y,z,t = q
    u = (F(0),F(1,3),F(2,3),F(1,3)+z,F(5,6)+z,F(2,3)+t,F(1,6)+t)
    v = (F(0),F(1,3)-x,F(2,3)-y,F(1,3),F(5,6)-x,F(2,3),F(1,6)-y)
    H = [[0]*7]
    for r in (u,v):
        assert all((modulus*a).denominator == 1 for a in r)
        H.append([int(modulus*a) % modulus for a in r])
    quotient = tuple((a-b) % modulus for a,b in zip(H[1],H[2]))
    expected = tuple(int(modulus*a) % modulus for a in (F(0),x,y,z,x+z,t,y+t))
    assert quotient == expected
    return H


def phi(q):
    x,y,z,t = q
    return tuple(-a % 1 for a in (z,t,x,y))


def c_decompositions(r, m):
    result = []
    for T in combinations(range(7),3):
        if {(r[j]-r[T[0]]) % m for j in T} != {0,m//3,2*m//3}:
            continue
        rest = set(range(7))-set(T)
        for P in combinations(sorted(rest),2):
            Q = tuple(sorted(rest-set(P)))
            if P > Q:
                continue
            if (r[P[0]]-r[P[1]]) % m == m//2 and (r[Q[0]]-r[Q[1]]) % m == m//2:
                result.append([list(T),list(P),list(Q)])
    return result


def b_decompositions(r,m):
    result=[]
    for T in combinations(range(7),5):
        P=tuple(j for j in range(7) if j not in T)
        if {(r[j]-r[T[0]]) % m for j in T} == {k*m//5 for k in range(5)} \
                and (r[P[0]]-r[P[1]]) % m == m//2:
            result.append([list(T),list(P)])
    return result


def dephase(H, rows, columns, m):
    b,r,s=rows
    a=columns[0]
    return tuple((H[k][j]-H[b][j]-H[k][a]+H[b][a]) % m
                 for k in (r,s) for j in columns)


def main():
    data=json.loads((ROOT/'enumeration.json').read_text())
    all_q = {}
    all_q['ACC'] = {tuple(F(a,7) for a in q) for q in product(range(7),repeat=4)
                    if len({0,q[0],q[1],q[2],(q[0]+q[2])%7,q[3],(q[1]+q[3])%7})==7}
    full_s=[(j,k) for j,k in product(range(1,5),repeat=2) if j!=k and (j+k)%5]
    all_q['BCC']=set()
    for j,k in full_s:
        for eps in (1,-1):
            h=(F(eps,4)-F(j+k,5)) % 1
            all_q['BCC'].update(((h,F(j,5),(h+F(1,2))%1,F(k,5)),
                                 (F(j,5),h,F(k,5),(h+F(1,2))%1)))
    reps = {'ACC': [(f'X{i+1}',tuple(F(a,7) for a in q)) for i,q in enumerate(ACC_REPS)],'BCC':[]}
    for placement in (1,2):
        for j,k in S:
            for eps in (1,-1):
                h=(F(eps,4)-F(j+k,5)) % 1
                q=(h,F(j,5),(h+F(1,2))%1,F(k,5)) if placement==1 else \
                    (F(j,5),h,F(k,5),(h+F(1,2))%1)
                reps['BCC'].append((f'Y{placement}_{j}{k}_{"plus" if eps==1 else "minus"}',q))
    output={'types':{},'visited_transformations':0}
    column_orders=list(permutations(range(7)))
    for kind,m in (('ACC',42),('BCC',60)):
        assert len(all_q[kind])==(24 if kind=='ACC' else 32)
        primary={tuple(tuple(F(x) for x in r) for r in H) for H in data['types'][kind]['points']}
        literal={tuple(tuple(F(a,m) for a in r) for r in M(q,m)[1:]) for q in all_q[kind]}
        assert literal==primary
        lookup={dephase(M(q,m),(0,1,2),tuple(range(7)),m):label for label,q in reps[kind]}
        assert len(lookup)==(12 if kind=='ACC' else 16)
        certificates=[]
        for q in sorted(all_q[kind]):
            H=M(q,m)
            C=[c_decompositions(r,m) for r in H[1:]]
            assert all(len(ds)==1 for ds in C)
            quotient=[(a-b)%m for a,b in zip(H[1],H[2])]
            if kind=='ACC':
                assert set(quotient)=={k*m//7 for k in range(7)}
                quotient_support=[list(range(7))]
            else:
                Q=b_decompositions(quotient,m)
                assert len(Q)==1
                assert Q[0][1] in ([1,3],[2,5])
                quotient_support=Q[0]
            assert phi(q) in all_q[kind] and phi(q)!=q
            # Exact universal equivalence: swap u,v and columns 1,4,6,2,5,3,7.
            assert dephase(H,(0,2,1),(0,3,5,1,4,2,6),m)== \
                dephase(M(phi(q),m),(0,1,2),tuple(range(7)),m)
            intended=[label for label,r in reps[kind] if q==r or phi(q)==r]
            assert len(intended)==1
            matches={}
            for rows in permutations(range(3)):
                for columns in column_orders:
                    output['visited_transformations']+=1
                    label=lookup.get(dephase(H,rows,columns,m))
                    if label is not None:
                        matches.setdefault(label,{'rows':[r+1 for r in rows],
                                                  'columns':[j+1 for j in columns]})
            assert set(matches)==set(intended),(kind,q,matches,intended)
            certificates.append({'parameter_phases':[str(a) for a in q],
                                 'class':intended[0],**matches[intended[0]],
                                 'u_support':[[j+1 for j in block] for block in C[0][0]],
                                 'v_support':[[j+1 for j in block] for block in C[1][0]],
                                 'quotient_support':[[j+1 for j in block] for block in quotient_support]})
        output['types'][kind]={'classes':len(reps[kind]),'points':len(all_q[kind]),
                              'matrix_modulus':m,
                              'representatives':[{'label':label,'parameter_phases':[str(a) for a in q],
                                                  'exponents':M(q,m)} for label,q in reps[kind]],
                              'certificates':certificates}
        print(kind,': exact formula coverage, unique C supports, involution and inequivalence verified.',flush=True)
    assert output['visited_transformations']==1693440
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print('Independent entrywise check:',output['visited_transformations'],'row/column transformations.',flush=True)


if __name__=='__main__':
    main()
