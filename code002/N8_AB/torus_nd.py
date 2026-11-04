"""Exact connected affine subtori of T^N, using integer Smith decomposition."""
N=3
from fractions import Fraction as F
from itertools import product
from math import lcm
import json
from pathlib import Path
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

def enc(x):
    if isinstance(x,F):return [x.numerator,x.denominator]
    if isinstance(x,(tuple,list)):return [enc(y) for y in x]
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    return x
def rat(x):return F(int(x.p),int(x.q)) if hasattr(x,'q') else F(x)
def torus(A,c):
    """A is saturated; c is the value of A*t. Canonical row basis."""
    if not len(A):return {'A':[],'c':[],'o':[F(0)]*N,'V':[[int(i==j) for j in range(N)] for i in range(N)],'dim':N}
    A=S.Matrix(A);c=S.Matrix([S.Rational(x.numerator,x.denominator) for x in c])
    B=hermite_normal_form(A.T).T
    # B = T*A, recover T from independent columns.
    inds=A.rref()[1];T=B[:,list(inds)]*A[:,list(inds)].inv();d=T*c
    dm=DomainMatrix.from_Matrix(B).convert_to(ZZ)
    sd,ud,vd=smith_normal_decomp(dm);D,U,V=[x.to_Matrix() for x in (sd,ud,vd)]
    r=B.rows;assert all(abs(D[i,i])==1 for i in range(r))
    y=S.zeros(N,1);z=U*d
    for i in range(r):y[i]=z[i]/D[i,i]
    o=V*y
    return {'A':[[int(x) for x in row] for row in B.tolist()], 'c':[rat(x)%1 for x in d],
            'o':[rat(x)%1 for x in o], 'V':[[int(x) for x in row] for row in V[:,r:].tolist()],'dim':N-r}
def key(t):return (tuple(map(tuple,t['A'])),tuple(t['c']))
def decode(t):
    return {**t,'c':[F(*x) for x in t['c']],'o':[F(*x) for x in t['o']]}
def split(eq,mod=30):
    if not eq:return [torus([],[])]
    A=S.Matrix([e[:N] for e in eq]);c=S.Matrix([S.Rational(-e[N],mod) for e in eq]);r=A.rank()
    dm=DomainMatrix.from_Matrix(A).convert_to(ZZ)
    sd,ud,vd=smith_normal_decomp(dm);D,U,V=[x.to_Matrix() for x in (sd,ud,vd)]
    assert U*A*V==D
    z=U*c
    if any(x.q!=1 for x in z[r:]):return []
    B=V.inv()[:r,:];out=[]
    for ns in product(*(range(abs(int(D[i,i]))) for i in range(r))):
        y=S.zeros(N,1)
        for i in range(r):y[i]=(z[i]+ns[i])/D[i,i]
        out.append(torus(B.tolist(),[rat(x) for x in (B*V*y)]))
    return out
def subset(t,u):
    for a,c in zip(u['A'],u['c']):
        if (sum(F(x)*y for x,y in zip(a,t['o']))-c)%1:return False
        for j in range(t['dim']):
            if sum(a[i]*t['V'][i][j] for i in range(N)):return False
    return True
def equations(t,mod):return [a+[int((-c*mod)%mod)] for a,c in zip(t['A'],t['c'])]
def intersection(t,u):
    m=lcm(30,*(x.denominator for x in t['c']+u['c']))
    return split(equations(t,m)+equations(u,m),m)
def specialize(H,t,mod):
    return [[[int((F(e[0],30)+sum(F(e[i+1])*t['o'][i] for i in range(N)))*mod)%mod]
             +[sum(e[i+1]*t['V'][i][j] for i in range(N)) for j in range(t['dim'])] for e in row] for row in H]
def lift(row,t,mod):
    if t['dim']:
        V=S.Matrix(t['V']);dm=DomainMatrix.from_Matrix(V).convert_to(ZZ)
        sd,ud,vd=smith_normal_decomp(dm);D,U,W=[x.to_Matrix() for x in (sd,ud,vd)]
        dd=t['dim'];J=S.zeros(dd,N)
        for i in range(dd):assert abs(D[i,i])==1;J[i,i]=1/D[i,i]
        L=W*J*U;assert L*V==S.eye(dd)
    out=[]
    for e in row:
        a=list(S.Matrix([e[1:]])*L) if t['dim'] else [0]*N
        a=[int(x) for x in a];c=(F(e[0],mod)-sum(x*y for x,y in zip(a,t['o'])))%1
        out.append([c]+a)
    return out
def write_scan(path,items):
    mod=lcm(30,*(x.denominator for it in items for x in it['torus']['o']))
    with Path(path).open('w') as f:
        f.write(str(len(items))+'\n')
        for idx,it in enumerate(items):
            t=it['torus'];f.write(f"{it['class_id']} {idx} {t['dim']}\n")
            for row in specialize(it['H'],t,mod):f.write(' '.join(str(x) for e in row for x in e)+'\n')
    return mod
