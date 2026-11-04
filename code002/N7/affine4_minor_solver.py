"""Independent rational-minor solver for integer congruences on a torus.

Does not use integer-unimodular elimination. Cached work depends only on A.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import deque
from math import lcm


def rref(A):
    A = [[F(x) for x in row] for row in A]
    pivots = []
    r = 0
    for j in range(len(A[0]) if A else 0):
        k = next((k for k in range(r, len(A)) if A[k][j]), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        d = A[r][j]
        A[r] = [x/d for x in A[r]]
        for k in range(len(A)):
            if k != r:
                d = A[k][j]
                A[k] = [x-d*y for x, y in zip(A[k], A[r])]
        pivots.append(j)
        r += 1
        if r == len(A):
            break
    return A, pivots


def det(A):
    if not A:
        return 1
    return sum((-1)**j*A[0][j]*det([r[:j]+r[j+1:] for r in A[1:]])
               for j in range(len(A)))


@lru_cache(None)
def prepare(A, n):
    _, columns = rref(A)
    rank = len(columns)
    if not rank:
        return (), (), (), ((F(0),)*0,), tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
    _, rows = rref(list(zip(*A)))
    rows = rows[:rank]
    M = [[A[i][j] for j in columns] for i in rows]
    aug = [r+[int(i == j) for j in range(rank)] for i, r in enumerate(M)]
    reduced, _ = rref(aug)
    inv = [r[rank:] for r in reduced]
    generators = [tuple(inv[i][j] % 1 for i in range(rank)) for j in range(rank)]
    zero = (F(0),)*rank
    group, todo = {zero}, deque([zero])
    while todo:
        x = todo.popleft()
        for g in generators:
            y = tuple((a+b) % 1 for a, b in zip(x, g))
            if y not in group:
                group.add(y)
                todo.append(y)
    assert len(group) == abs(det(M))
    free = [j for j in range(n) if j not in columns]
    B = [[F(0)]*len(free) for _ in range(n)]
    for k, j in enumerate(free):
        B[j][k] = F(1)
        for h, p in enumerate(columns):
            B[p][k] = -sum(inv[h][s]*A[i][j] for s, i in enumerate(rows))
    for k in range(len(free)):
        scale = lcm(*(B[i][k].denominator for i in range(n)))
        for i in range(n):
            B[i][k] = int(scale*B[i][k])
    B = tuple(tuple(r) for r in B)
    assert all(sum(ar[i]*B[i][j] for i in range(n)) == 0
               for ar in A for j in range(len(free)))
    return tuple(rows), tuple(columns), inv, tuple(sorted(group)), B


def solve(A, b, n=4):
    # Reject incompatible repeated equations, including the zero equation.
    known = {}
    for ar, target in zip(A, b):
        ar = tuple(ar)
        if not any(ar):
            if target % 1:
                return []
            continue
        sign = 1 if next(x for x in ar if x) > 0 else -1
        key = tuple(sign*x for x in ar)
        value = (sign*target) % 1
        if key in known and known[key] != value:
            return []
        known[key] = value
    A = tuple(tuple(r) for r in A)
    rows, columns, inv, group, B = prepare(A, n)
    origin = [sum(inv[h][s]*b[i] for s, i in enumerate(rows)) for h in range(len(columns))]
    result = []
    for shift in group:
        c = [F(0)]*n
        for h, j in enumerate(columns):
            c[j] = (origin[h]+shift[h]) % 1
        if all((sum(a*x for a, x in zip(ar, c))-target) % 1 == 0 for ar, target in zip(A, b)):
            result.append((tuple(c), B))
    return result
