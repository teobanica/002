"""Exact phase and Smith utilities for symbolic completion families."""
from pathlib import Path
import json,sys
from fractions import Fraction as F
from itertools import combinations,product
from collections import Counter
import sympy as S
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from exact_rows import profiles,quotient

R=Path(__file__).resolve().parents[1]
def read(n):return json.loads((R/'outputs'/n).read_text())
def decode(row):return tuple((F(e[0],30),)+tuple(map(F,e[1:])) for e in row)
def pad(row,n):return tuple(e+(F(0),)*n for e in row)
def enc(o):
    if isinstance(o,F):return [o.numerator,o.denominator]
    if isinstance(o,(list,tuple,set)):return [enc(x) for x in sorted(o)] if isinstance(o,set) else [enc(x) for x in o]
    if isinstance(o,dict):return {k:enc(v) for k,v in o.items()}
    return o
def substitute(row,base,values,extra=0):
    result=[]
    for e in row:
        v=list(e[:base+1]) + [F(0)]*extra
        for k,value in enumerate(values):
            for j,x in enumerate(value):v[j]+=e[base+1+k]*x
        v[0]%=1;result.append(tuple(v))
    return tuple(result)
def parameter_for(row,target,base):
    j=next(j for j,e in enumerate(row) if e[base+1])
    k=row[j][base+1]
    assert k==1
    value=tuple((target[j][i]-row[j][i]) if i else (target[j][0]-row[j][0])%1 for i in range(base+1))
    return value if substitute(row,base,[value])==target else None
def solutions(eq,f,base):
    if not eq:
        return [[(F(0),)*(base+1)+tuple(F(i==j) for j in range(f)) for i in range(f)]]
    A=S.Matrix([e[:f] for e in eq]);B=S.Matrix([e[f:f+base]+[S.Rational(e[-1],30)] for e in eq])
    sd,ud,vd=smith_normal_decomp(DomainMatrix.from_Matrix(A).convert_to(ZZ))
    D,U,V=[x.to_Matrix() for x in (sd,ud,vd)];assert U*A*V==D
    rank=A.rank();C=-U*B
    assert all(x==0 for x in C[rank:,:])
    out=[]
    for ks in product(*(range(abs(int(D[i,i]))) for i in range(rank))):
        Y=S.zeros(f,base+1+f-rank)
        for i in range(rank):
            for j in range(base+1):Y[i,j]=(C[i,j]+(ks[i] if j==base else 0))/D[i,i]
        for i in range(rank,f):Y[i,base+1+i-rank]=1
        X=V*Y
        out.append([tuple([F(X[i,base])%1]+[F(X[i,j]) for j in range(base)]+[F(X[i,j]) for j in range(base+1,X.cols)]) for i in range(f)])
    return out
def cliques(adj,need):
    out=[]
    def rec(xs,sol):
        if len(sol)==need:out.append(sol);return
        while len(xs)+len(sol)>=need:
            i=xs[0];xs=xs[1:];rec([j for j in xs if adj[i][j]],sol+[i])
    rec(list(range(len(adj))),[]);return out

