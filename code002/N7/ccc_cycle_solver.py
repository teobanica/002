"""Exact 3+2+2 cycle-coupling enumeration; Python standard library only.

Angles are measured in turns, in Q/Z; free parameters range over R/Z.
No floating point arithmetic and no bound on parameter root orders is used.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import gcd
from collections import Counter
import json
from pathlib import Path
from time import perf_counter

from support_orbits import support_partitions, BASE


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


def cycle_row(partition, orientation, variable_offset):
    """Row with ternary anchor 1 and independent phases on the two pairs."""
    c = [F(0)]*7
    B = [[0]*4 for _ in range(7)]
    triple, pair_a, pair_b = partition
    for j, col in enumerate(triple):
        c[col] = F(orientation*j, 3) % 1
    for j, pair in enumerate((pair_a, pair_b)):
        for col in pair:
            B[col][variable_offset+j] = 1
        c[pair[1]] = F(1, 2)
    return c, B


def quotient_constraints(u, v, partition, orientation):
    uc, UB = u
    vc, VB = v
    c = [x-y for x, y in zip(uc, vc)]
    B = [[x-y for x, y in zip(a, b)] for a, b in zip(UB, VB)]
    triple, pair_a, pair_b = partition
    edges = [(triple[0], triple[1], F(orientation, 3)),
             (triple[0], triple[2], F(2*orientation, 3)),
             (pair_a[0], pair_a[1], F(1, 2)),
             (pair_b[0], pair_b[1], F(1, 2))]
    A = [tuple(y-x for x, y in zip(B[i], B[j])) for i, j, _ in edges]
    rhs = [(target-c[j]+c[i]) % 1 for i, j, target in edges]
    return A, rhs


def substitute_row(row, solution):
    c, B = row
    x, X = solution
    d = len(X[0])
    constant = [(ci+dot(bi, x)) % 1 for ci, bi in zip(c, B)]
    basis = [[sum(bi[k]*X[k][j] for k in range(4)) for j in range(d)] for bi in B]
    # Normalize entry 0 of each row, retaining the all-ones extra row.
    return (tuple((ci-constant[0]) % 1 for ci in constant),
            tuple(tuple(x-y for x, y in zip(bi, basis[0])) for bi in basis))


def verify_solution(A, b, c, B):
    for a, target in zip(A, b):
        assert (dot(a, c)-target) % 1 == 0
        assert all(sum(a[i]*B[i][j] for i in range(len(a))) == 0
                   for j in range(len(B[0])))


def enumerate_components():
    parts = support_partitions()
    oriented = list(product(range(len(parts)), (1, -1)))
    u = cycle_row(BASE, 1, 0)
    counts = Counter()
    records = []
    for vi, ve in oriented:
        v = cycle_row(parts[vi], ve, 2)
        for zi, ze in oriented:
            A, rhs = quotient_constraints(u, v, parts[zi], ze)
            solutions = solve_torus(A, rhs, 4)
            counts['systems'] += 1
            counts['consistent_systems'] += bool(solutions)
            for solution in solutions:
                c, B = solution
                verify_solution(A, rhs, c, B)
                d = len(B[0])
                counts[f'components_dim_{d}'] += 1
                rows = [substitute_row(row, solution) for row in (u, v)]
                records.append({'support': [vi, ve, zi, ze], 'dimension': d,
                                'constant': [[str(x) for x in r[0]] for r in rows],
                                'basis': [r[1] for r in rows]})
    return dict(counts), records


def self_check():
    # 2x=0 has two torsion components; incompatible equations have none.
    assert {c for c, _ in solve_torus([[2]], [F(0)], 1)} == {(F(0),), (F(1, 2),)}
    assert not solve_torus([[2], [2]], [F(0), F(1, 2)], 1)
    for r in product(range(-2, 3), repeat=4):
        if not any(r):
            continue
        V, g = unimodular_reduce(r)
        assert [sum(r[i]*V[i][j] for i in range(4)) for j in range(4)] == [g, 0, 0, 0]
        assert g == reduce(gcd, r)


if __name__ == '__main__':
    started = perf_counter()
    self_check()
    counts, records = enumerate_components()
    output = Path(__file__).with_name('cycle_components.json')
    output.write_text(json.dumps({'counts': counts, 'components': records}, indent=2)+'\n')
    print(json.dumps(counts, indent=2))
    print(f'Elapsed: {perf_counter()-started:.2f} seconds; saved {len(records)} components.')
