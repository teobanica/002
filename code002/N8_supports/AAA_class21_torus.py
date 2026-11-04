from fractions import Fraction as F
from math import gcd,lcm
from itertools import combinations
def egcd(a,b):
    if not b:return (abs(a),1 if a>=0 else -1,0)
    g,u,v=egcd(b,a%b);return g,v,u-(a//b)*v
def solve(eq):
    """Equations A*x+B*y=c mod 1. Return full, primitive lines, or points."""
    eq=[(int(a),int(b),F(c)%1) for a,b,c in eq]
    if any(not(a or b) and c for a,b,c in eq):return False,[],[]
    nz=[e for e in eq if e[0] or e[1]]
    if not nz:return True,[],[]
    a,b,c=nz[0];g=gcd(a,b);a//=g;b//=g
    _,u,v=egcd(a,b);lines=[];points=[]
    for k in range(g):
        d=((c+k)/g)%1;x=u*d;y=v*d
        if any(not(A*b-B*a) and (A*x+B*y-C)%1 for A,B,C in nz):continue
        pp=next(((A,B,C) for A,B,C in nz if A*b-B*a),None)
        if pp is None:
            aa,bb=a,b;dd=d
            if aa<0 or (aa==0 and bb<0):aa=-aa;bb=-bb;dd=(-dd)%1
            lines.append((aa,bb,dd));continue
        A,B,C=pp;h=A*b-B*a
        for k2 in range(abs(h)):
            t=(C-A*x-B*y+k2)/h;xx=(x+b*t)%1;yy=(y-a*t)%1
            if all((A*xx+B*yy-C)%1==0 for A,B,C in nz):points.append((xx,yy))
    return False,sorted(set(lines)),sorted(set(points))
def matchings(vertices=tuple(range(8))):
    if not vertices:yield ();return
    a=vertices[0]
    for b in vertices[1:]:
        for rest in matchings(tuple(j for j in vertices[1:] if j!=b)):yield ((a,b),)+rest
MATCH=list(matchings())
def line_param(line):
    a,b,c=line;g,u,v=egcd(a,b);assert g==1
    return (u*c,v*c),(b,-a)
def on_line(p,l):return (l[0]*p[0]+l[1]*p[1]-l[2])%1==0
def encode(o):
    if isinstance(o,F):return [o.numerator,o.denominator]
    if isinstance(o,(list,tuple)):return [encode(x) for x in o]
    if isinstance(o,dict):return {k:encode(v) for k,v in o.items()}
    return o
def decode_line(l):return l[0],l[1],F(*l[2])
def decode_point(p):return F(*p[0]),F(*p[1])
