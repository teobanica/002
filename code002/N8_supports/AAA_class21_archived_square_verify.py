#!/usr/bin/env python3
"""Exact symbolic verification of the two 5-parameter A-only 8x8 families.
Each entry is represented by (exponent vector in x,y,z,w,t, sign parity).
No floating point is used.
"""
from itertools import combinations
from collections import Counter

def M(x=0,y=0,z=0,w=0,t=0,h=0): return ((x,y,z,w,t),h&1)
def diff(r,s): return tuple((tuple(a-b for a,b in zip(c1,c2)),(h1-h2)&1) for (c1,h1),(c2,h2) in zip(r,s))
def is_A(q):
    C=Counter(q); coeff={c for c,h in q}
    return all(C[(c,0)]==C[(c,1)] for c in coeff)
def check(H):
    assert len(H)==8 and all(len(r)==8 for r in H)
    assert all(is_A(diff(H[i],H[j])) for i,j in combinations(range(8),2))
    HT=list(zip(*H))
    assert all(is_A(diff(HT[i],HT[j])) for i,j in combinations(range(8),2))

H0=[
 [M() for _ in range(8)],
 [M(),M(h=1),M(x=1),M(x=1,h=1),M(x=1),M(x=1,h=1),M(),M(h=1)],
 [M(),M(y=1,h=1),M(z=1),M(w=1,h=1),M(z=1,h=1),M(w=1),M(h=1),M(y=1)],
 [M(),M(t=1),M(y=-1,z=1,t=1),M(w=1),M(y=-1,z=1,t=1,h=1),M(w=1,h=1),M(h=1),M(t=1,h=1)],
 [M(),M(t=1,h=1),M(y=-1,z=1,t=1,h=1),M(w=1),M(y=-1,z=1,t=1),M(w=1,h=1),M(h=1),M(t=1)],
 [M(),M(y=1),M(z=1,h=1),M(w=1,h=1),M(z=1),M(w=1),M(h=1),M(y=1,h=1)],
 [M(),M(),M(h=1),M(h=1),M(h=1),M(h=1),M(),M()],
 [M(),M(h=1),M(x=1,h=1),M(x=1),M(x=1,h=1),M(x=1),M(),M(h=1)],
]
H1=[
 [M() for _ in range(8)],
 [M(),M(h=1),M(x=1,h=1),M(x=1),M(y=1),M(y=1,h=1),M(),M(h=1)],
 [M(),M(z=1,h=1),M(w=1,h=1),M(w=1,h=1),M(w=1),M(w=1),M(h=1),M(z=1)],
 [M(),M(z=1),M(t=1),M(t=1,h=1),M(x=-1,y=1,t=1),M(x=-1,y=1,t=1,h=1),M(h=1),M(z=1,h=1)],
 [M(),M(z=1),M(t=1,h=1),M(t=1),M(x=-1,y=1,t=1,h=1),M(x=-1,y=1,t=1),M(h=1),M(z=1,h=1)],
 [M(),M(z=1,h=1),M(w=1),M(w=1),M(w=1,h=1),M(w=1,h=1),M(h=1),M(z=1)],
 [M(),M(),M(h=1),M(h=1),M(h=1),M(h=1),M(),M()],
 [M(),M(h=1),M(x=1),M(x=1,h=1),M(y=1,h=1),M(y=1),M(),M(h=1)],
]
check(H0); check(H1)
print('verified H0 and H1: all 28 row quotients and all 28 column quotients have profile A')
