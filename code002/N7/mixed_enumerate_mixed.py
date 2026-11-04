"""Exhaust all ACC/BCC supports with unrestricted real phases modulo one.

Every congruence system is solved independently in two ways. Neither solver
assumes a root order. A second pass can use all 105 supports without symmetry.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter
from pathlib import Path
from math import lcm
import json
from torus_solver import solve_torus
from minor_solver import solve as minor_solve

ROOT = Path(__file__).resolve().parent
BASE = ((0, 1, 2), (3, 4), (5, 6))


def partitions():
    out = []
    for T in combinations(range(7), 3):
        a, b, c, d = [j for j in range(7) if j not in T]
        for P, Q in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
            out.append((T, P, Q))
    assert len(set(out)) == 105
    return out


def supports():
    group = []
    for shift, swap, f, g in product(range(3), range(2), range(2), range(2)):
        T = tuple((j+shift) % 3 for j in range(3))
        P = (4, 3) if f else (3, 4)
        Q = (6, 5) if g else (5, 6)
        group.append(T+(Q+P if swap else P+Q))
    assert len(set(group)) == 24
    remaining = set(partitions())
    answer = []
    while remaining:
        P = min(remaining)
        orbit = set()
        for g in group:
            T = tuple(sorted(g[j] for j in P[0]))
            pairs = sorted(tuple(sorted(g[j] for j in pair)) for pair in P[1:])
            orbit.add((T, *pairs))
        assert orbit <= remaining
        answer.append((P, len(orbit)))
        remaining -= orbit
    assert len(answer) == 9 and sum(n for _, n in answer) == 105
    return answer


def row(P, sign, offset):
    c = [F(0)]*7
    B = [[0]*4 for _ in range(7)]
    for k, j in enumerate(P[0]):
        c[j] = F(sign*k, 3) % 1
    for k, pair in enumerate(P[1:]):
        for j in pair:
            B[j][offset+k] = 1
        c[pair[1]] = F(1, 2)
    return c, B


def patterns(kind):
    if kind == 'ACC':
        return [tuple((0, j, F(k, 7)) for j, k in zip(range(1, 7), p))
                for p in permutations(range(1, 7))]
    out = []
    for T in combinations(range(7), 5):
        pair = tuple(j for j in range(7) if j not in T)
        for p in permutations(range(1, 5)):
            out.append(tuple((T[0], j, F(k, 5)) for j, k in zip(T[1:], p))
                       +((*pair, F(1, 2)),))
    assert len(out) == 504
    return out


def equations(U, V, edges):
    c = [a-b for a, b in zip(U[0], V[0])]
    B = [[a-b for a, b in zip(r, s)] for r, s in zip(U[1], V[1])]
    A = [[B[j][k]-B[i][k] for k in range(4)] for i, j, _ in edges]
    b = [(t-c[j]+c[i]) % 1 for i, j, t in edges]
    return A, b


def evaluate(U, V, point):
    result = []
    for c, B in (U, V):
        r = [(a+sum(x*y for x, y in zip(b, point))) % 1 for a, b in zip(c, B)]
        result.append(tuple((x-r[0]) % 1 for x in r))
    return tuple(result)


def canonical(H):
    H = ((F(0),)*7,)+tuple(H)
    choices = []
    for b, r, s in permutations(range(3)):
        pairs = [((H[r][j]-H[b][j]) % 1, (H[s][j]-H[b][j]) % 1) for j in range(7)]
        for x0, y0 in pairs:
            choices.append(tuple(sorted(((x-x0) % 1, (y-y0) % 1) for x, y in pairs)))
    return min(choices)


def main():
    U = row(BASE, 1, 0)
    report = {'supports': [{'case': i+1, 'partition': P, 'orbit_size': n}
                           for i, (P, n) in enumerate(supports())], 'types': {}}
    for kind in ('ACC', 'BCC'):
        by_case, records, points = [], [], set()
        for case, (VP, orbit_size) in enumerate(supports(), 1):
            count = Counter()
            for sign in (1, -1):
                V = row(VP, sign, 2)
                for pattern, edges in enumerate(patterns(kind)):
                    A, b = equations(U, V, edges)
                    first = solve_torus(A, b, 4)
                    second = minor_solve(A, b, 4)
                    # If any system has a continuous solution, stop rather
                    # than silently forcing the planned finite classification.
                    assert all(len(B[0]) == 0 for _, B in first+second), (kind, case, sign, pattern)
                    X = {tuple(x % 1 for x in c) for c, _ in first}
                    Y = {tuple(x % 1 for x in c) for c, _ in second}
                    assert X == Y, (kind, case, sign, pattern, X, Y)
                    count['systems'] += 1
                    count['consistent_systems'] += bool(X)
                    count['solutions'] += len(X)
                    for x in sorted(X):
                        assert all((sum(a*y for a, y in zip(ar, x))-t) % 1 == 0
                                   for ar, t in zip(A, b))
                        assert case == 8 and sign == 1
                        H = evaluate(U, V, x)
                        points.add(H)
                        records.append({'case': case, 'orientation': sign, 'pattern': pattern,
                                        'phase_solution': [str(a) for a in x],
                                        'rows': [[str(a) for a in r] for r in H]})
            by_case.append({'case': case, 'orbit_size': orbit_size, **count})
            print(kind, 'case', case, dict(count), flush=True)
        modulus = lcm(*(x.denominator for H in points for r in H for x in r))
        quotient_modulus = lcm(*((a-b).denominator for H in points for a, b in zip(*H)))
        keys = {canonical(H) for H in points}
        assert len(keys) == (12 if kind == 'ACC' else 16)
        assert modulus == (42 if kind == 'ACC' else 60)
        assert quotient_modulus == (7 if kind == 'ACC' else 20)
        report['types'][kind] = {'by_case': by_case, 'records': records,
                                'point_count': len(points), 'classes': len(keys),
                                'matrix_modulus': modulus, 'quotient_modulus': quotient_modulus,
                                'points': [[[str(a) for a in r] for r in H] for H in sorted(points)]}
        print(kind, 'distinct points', len(points), 'classes', len(keys), 'root order', modulus, flush=True)
    (ROOT/'enumeration.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Both exact solvers agree on all 22032 systems.', flush=True)


if __name__ == '__main__':
    main()
