"""Exact integer-unimodular elimination on R/Z; no root-order bound."""
from fractions import Fraction as F


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def egcd(a, b):
    """Return g >= 0 and s,t with s*a+t*b=g."""
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r-q*r
        old_s, s = s, old_s-q*s
        old_t, t = t, old_t-q*t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def unimodular_reduce(r):
    """Return V, g with r V=(g,0,...,0), det(V)=+-1, g=gcd(r)."""
    r = list(r)
    n = len(r)
    V = [[int(i == j) for j in range(n)] for i in range(n)]
    for j in range(1, n):
        a, b = r[0], r[j]
        if b == 0:
            continue
        g, s, t = egcd(a, b)
        for row in V:
            x, y = row[0], row[j]
            row[0], row[j] = s*x+t*y, -(b//g)*x+(a//g)*y
        r[0], r[j] = g, 0
    if r[0] < 0:
        r[0] = -r[0]
        for row in V:
            row[0] = -row[0]
    return V, r[0]


def solve_torus(A, b, n):
    """All connected components of A*x=b mod 1 as x=c+B*t.

    A is integer, b rational. B is integer; t ranges freely over a torus.
    Integer unimodular variable changes are bijections on the torus.
    The branches (rhs+k)/g, 0<=k<g, exhaust g*t=rhs mod 1.
    """
    states = [(tuple(F(0) for _ in range(n)),
               tuple(tuple(int(i == j) for j in range(n)) for i in range(n)))]
    for a, target in zip(A, b):
        following = []
        for c, B in states:
            d = len(B[0])
            r = [sum(a[i]*B[i][j] for i in range(n)) for j in range(d)]
            rhs = (target-dot(a, c)) % 1
            if not any(r):
                if rhs == 0:
                    following.append((c, B))
                continue
            V, g = unimodular_reduce(r)
            BV = [[sum(B[i][k]*V[k][j] for k in range(d))
                   for j in range(d)] for i in range(n)]
            new_B = tuple(tuple(row[1:]) for row in BV)
            for k in range(g):
                new_c = tuple((c[i]+BV[i][0]*(rhs+k)/g) % 1 for i in range(n))
                following.append((new_c, new_B))
        states = following
        if not states:
            break
    return states
