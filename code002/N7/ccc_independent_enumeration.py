"""Independent rational solver and fresh CCC support enumeration.

No floating point arithmetic. The second solver uses a nonsingular minor,
rational inversion and a finite lattice quotient, not torus elimination.
"""
from fractions import Fraction as F
from itertools import combinations, product, permutations
from math import gcd, lcm
from collections import Counter, deque
from pathlib import Path
import json, time

ROOT=Path(__file__).parent
ZERO=F(0)
BASE=((0,1,2),(3,4),(5,6))

def partitions():
    out=[]
    for T in combinations(range(7),3):
        a,b,c,d=[j for j in range(7) if j not in T]
        for P,Q in [((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))]:
            out.append((T,P,Q))
    assert len(out)==len(set(out))==105
    return out

def stabilizer():
    out=[]
    for shift,swap,flip1,flip2 in product(range(3),range(2),range(2),range(2)):
        T=tuple((i+shift)%3 for i in range(3))
        P=(4,3) if flip1 else (3,4)
        Q=(6,5) if flip2 else (5,6)
        out.append(T+(Q+P if swap else P+Q))
    assert len(set(out))==24
    return out

def act(P,g):
    T,A,B=P
    return (tuple(sorted(g[i] for i in T)),
            *sorted([tuple(sorted(g[i] for i in A)),tuple(sorted(g[i] for i in B))]))

def representatives():
    left=set(partitions()); out=[]
    for P in sorted(left.copy()):
        if P not in left: continue
        O={act(P,g) for g in stabilizer()}
        assert O<=left
        out.append((P,len(O))); left-=O
    assert not left and len(out)==9
    assert sum(n for _,n in out)==105
    return out

def row(P,sgn,which):
    # Affine phases in four independent variables a,b,c,d.
    const=[ZERO]*7; coeff=[[0]*4 for _ in range(7)]
    for k,j in enumerate(P[0]): const[j]=F(sgn*k,3)%1
    for k,pair in enumerate(P[1:]):
        for j in pair: coeff[j][2*which+k]=1
        const[pair[1]]=F(1,2)
    return const,coeff

def equations(U,V,P,sgn):
    c=[x-y for x,y in zip(U[0],V[0])]
    B=[[x-y for x,y in zip(a,b)] for a,b in zip(U[1],V[1])]
    T,P1,P2=P
    constraints=[(T[0],T[1],F(sgn,3)),(T[0],T[2],F(2*sgn,3)),
                 (*P1,F(1,2)),(*P2,F(1,2))]
    A=[]; rhs=[]
    for i,j,value in constraints:
        A.append([B[j][k]-B[i][k] for k in range(4)])
        rhs.append((value-c[j]+c[i])%1)
    return A,rhs

def rref(A):
    R=[[F(x) for x in row] for row in A]
    if not R: return R,[]
    piv=[]; p=0
    for j in range(len(R[0])):
        k=next((k for k in range(p,len(R)) if R[k][j]),None)
        if k is None: continue
        R[p],R[k]=R[k],R[p]
        a=R[p][j]; R[p]=[v/a for v in R[p]]
        for k in range(len(R)):
            if k!=p and R[k][j]:
                a=R[k][j]; R[k]=[x-a*y for x,y in zip(R[k],R[p])]
        piv.append(j); p+=1
        if p==len(R): break
    return R,piv

def inverse(A):
    n=len(A)
    R,p=rref([list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(A)])
    assert p[:n]==list(range(n))
    return [row[n:] for row in R]

def determinant(A):
    n=len(A)
    if n==0:return 1
    if n==1:return A[0][0]
    return sum((-1)**j*A[0][j]*determinant([r[:j]+r[j+1:] for r in A[1:]]) for j in range(n))

def matmul(A,B):
    return [[sum(x*y for x,y in zip(a,b)) for b in zip(*B)] for a in A]

def minor_solve(A,b,n=4):
    """Exhaust the solutions by a nonsingular rank-sized minor.
    Integer basis columns may parametrize a torus noninjectively.
    """
    _,columns=rref(A)
    rank=len(columns)
    if rank==0:
        return [] if any(x%1 for x in b) else [([ZERO]*n,[[int(i==j) for j in range(n)] for i in range(n)])]
    _,rows=rref(list(zip(*A)))
    rows=rows[:rank]
    M=[[A[i][j] for j in columns] for i in rows]
    inv=inverse(M); D=abs(determinant(M))
    generators=[tuple(inv[i][j]%1 for i in range(rank)) for j in range(rank)]
    start=(ZERO,)*rank
    group={start}; queue=deque([start])
    while queue:
        x=queue.popleft()
        for g in generators:
            y=tuple((a+c)%1 for a,c in zip(x,g))
            if y not in group: group.add(y);queue.append(y)
    assert len(group)==D
    free=[j for j in range(n) if j not in columns]
    rational_basis=[[ZERO]*len(free) for _ in range(n)]
    for k,j in enumerate(free):
        rational_basis[j][k]=F(1)
        for h,p in enumerate(columns):
            rational_basis[p][k]=-sum(inv[h][s]*A[i][j] for s,i in enumerate(rows))
    B=[[0]*len(free) for _ in range(n)]
    for k in range(len(free)):
        scale=lcm(*(rational_basis[i][k].denominator for i in range(n)))
        for i in range(n): B[i][k]=int(scale*rational_basis[i][k])
    origin=[sum(inv[h][s]*b[i] for s,i in enumerate(rows)) for h in range(rank)]
    solutions=[]
    for shift in sorted(group):
        c=[ZERO]*n
        for h,j in enumerate(columns): c[j]=(origin[h]+shift[h])%1
        if all((sum(a*x for a,x in zip(ar,c))-target)%1==0 for ar,target in zip(A,b)):
            assert all(sum(ar[i]*B[i][j] for i in range(n))==0 for ar in A for j in range(len(free)))
            solutions.append((c,B))
    return solutions

def substitute(U,V,solution):
    c,B=solution; d=len(B[0])
    constants=[]; basis=[]
    for rc,R in (U,V):
        C=[(rc[j]+sum(R[j][k]*c[k] for k in range(4)))%1 for j in range(7)]
        S=[[sum(R[j][k]*B[k][h] for k in range(4)) for h in range(d)] for j in range(7)]
        constants.extend((x-C[0])%1 for x in C)
        basis.extend(tuple(x-y for x,y in zip(s,S[0])) for s in S)
    return constants,basis

def graph_canonical(c,B):
    """Canonical affine torus with an injective coordinate chart.
    Also accepts nonprimitive bases from minor_solve; integrality is checked.
    """
    d=len(B[0])
    if d==0:return tuple(x%1 for x in c),tuple(() for _ in c),()
    for P in combinations(range(len(c)),d):
        M=[list(B[i]) for i in P]
        if not determinant(M):continue
        C=matmul(B,inverse(M))
        if any(x.denominator!=1 for r in C for x in r):continue
        C=tuple(tuple(int(x) for x in r) for r in C)
        residual=tuple((c[i]-sum(C[i][j]*c[p] for j,p in enumerate(P)))%1 for i in range(len(c)))
        return residual,C,P
    raise AssertionError('No integral coordinate chart; generalize canonicalization.')

def encode(f):
    c,B,P=f
    return {'constant':[str(x) for x in c],'basis':B,'pivot':P}

def decode(f):
    return tuple(F(x) for x in f['constant']),tuple(tuple(r) for r in f['basis']),tuple(f['pivot'])

def enumerate_all():
    from cycle_solver import solve_torus
    U=row(BASE,1,0); parts=partitions()
    records=[]; unique=set(); totals=Counter(); bycase=[]
    for case,(VP,size) in enumerate(representatives(),1):
        counts=Counter()
        for vs,zi,zs in product((1,-1),range(105),(1,-1)):
            V=row(VP,vs,1)
            A,b=equations(U,V,parts[zi],zs)
            primary=solve_torus(A,b,4)
            secondary=minor_solve(A,b)
            F1={graph_canonical(*substitute(U,V,s)) for s in primary}
            F2={graph_canonical(*substitute(U,V,s)) for s in secondary}
            assert F1==F2,(case,vs,zi,zs)
            counts['systems']+=1
            counts['consistent_systems']+=bool(F1)
            counts['raw_connected_components']+=len(primary)
            for f in F1:
                counts[f'dimension_{len(f[2])}']+=1
                unique.add(f)
            records.append({'case':case,'v_sign':vs,'z_partition':zi,'z_sign':zs,
                            'solutions':[encode(f) for f in sorted(F1)]})
        totals.update(counts)
        bycase.append({'case':case,'support':VP,'orbit_size':size,**counts})
        print('case',case,dict(counts),flush=True)
    result={'counts':dict(totals),'support_cases':bycase,
            'distinct_normalized_components':len(unique),
            'components':[encode(f) for f in sorted(unique)],'systems':records}
    (ROOT/'enumeration.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print('TOTAL',dict(totals),'distinct',len(unique),flush=True)
    return result

if __name__=='__main__':
    t=time.perf_counter()
    enumerate_all()
    print('Two independent solvers agree on every system;',round(time.perf_counter()-t,2),'seconds.')

