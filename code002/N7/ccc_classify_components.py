"""Exact containment and monomial-equivalence checks for cycle_solver.py.

Constants use Z/24Z, purely as an exact representation of the constants
returned by the solver; free parameters remain arbitrary circle variables.
"""
import json
from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter
from pathlib import Path
from time import perf_counter
from support_orbits import support_partitions

MOD = 24


def det(M):
    if not M:
        return 1
    if len(M) == 1:
        return M[0][0]
    return sum((-1)**j*M[0][j]*det([r[:j]+r[j+1:] for r in M[1:]])
               for j in range(len(M)))


def inv_unimodular(M):
    d = len(M)
    determinant = det(M)
    assert abs(determinant) == 1
    return [[(-1)**(i+j)*det([row[:i]+row[i+1:] for k, row in enumerate(M) if k != j])
             // determinant for j in range(d)] for i in range(d)]


def canonical(c, B):
    """Unique family representation with an intrinsic unimodular pivot.

    All encountered connected tori have a coordinate projection of degree 1;
    assert this property rather than assuming it for arbitrary input.
    """
    d = len(B[0])
    if not d:
        return (tuple(x % MOD for x in c), tuple(() for _ in c), ())
    pivot = next((p for p in combinations(range(len(B)), d)
                  if abs(det([list(B[i]) for i in p])) == 1), None)
    assert pivot is not None, 'No unimodular coordinate chart; extend algorithm.'
    inverse = inv_unimodular([list(B[i]) for i in pivot])
    C = tuple(tuple(sum(row[k]*inverse[k][j] for k in range(d)) for j in range(d))
              for row in B)
    residual = tuple((x-sum(row[k]*c[pivot[k]] for k in range(d))) % MOD
                     for x, row in zip(c, C))
    return residual, C, pivot


def dimension(family):
    return len(family[2])


def contains(big, small):
    bc, B, pivot = big
    sc, S, _ = small
    if len(S[0]) > len(pivot):
        return False
    for i, row in enumerate(B):
        if (sc[i]-bc[i]-sum(v*sc[p] for v, p in zip(row, pivot))) % MOD:
            return False
        for j, value in enumerate(S[i]):
            if value != sum(v*S[p][j] for v, p in zip(row, pivot)):
                return False
    return True


def from_record(record):
    c = [F(x)*MOD for row in record['constant'] for x in row]
    assert all(x.denominator == 1 for x in c)
    B = [tuple(row) for rr in record['basis'] for row in rr]
    return canonical([int(x) for x in c], B)


def oriented_cycles(c, B):
    """Every cycle support valid identically on this symbolic row."""
    answer = []
    for T, P, Q in support_partitions():
        if B[P[0]] != B[P[1]] or (c[P[1]]-c[P[0]]) % MOD != MOD//2:
            continue
        if B[Q[0]] != B[Q[1]] or (c[Q[1]]-c[Q[0]]) % MOD != MOD//2:
            continue
        if not (B[T[0]] == B[T[1]] == B[T[2]]):
            continue
        for sign in (1, -1):
            if all((c[T[k]]-c[T[0]]) % MOD == (sign*k*MOD//3) % MOD for k in (1, 2)):
                answer.append(((T, P, Q), sign))
    return answer


def normalize_rows(c, B, col_perm):
    return canonical(
        [(c[7*r+j]-c[7*r+col_perm[0]]) % MOD for r in range(2) for j in col_perm],
        [tuple(x-y for x, y in zip(B[7*r+j], B[7*r+col_perm[0]]))
         for r in range(2) for j in col_perm])


def normalized_orbit(family):
    """All equivalent families with first row (1,w,w²,a,-a,b,-b).

    Append the all-ones row; permute its three rows; dephase columns.
    Then use each identically valid ternary/pair decomposition of new row 1.
    This exhausts equivalences between normalized GENERIC families.
    """
    c, B, _ = family
    d = len(B[0])
    rows_c = [(0,)*7, c[:7], c[7:]]
    rows_B = [((0,)*d,)*7, B[:7], B[7:]]
    result = set()
    for base, first, second in permutations(range(3)):
        cc = tuple((rows_c[r][j]-rows_c[base][j]) % MOD
                   for r in (first, second) for j in range(7))
        BB = tuple(tuple(x-y for x, y in zip(rows_B[r][j], rows_B[base][j]))
                   for r in (first, second) for j in range(7))
        for (T, P, Q), sign in oriented_cycles(cc[:7], BB[:7]):
            oriented_T = T if sign == 1 else (T[0], T[2], T[1])
            for shift, swap, revP, revQ in product(range(3), range(2), range(2), range(2)):
                TT = oriented_T[shift:]+oriented_T[:shift]
                PP, QQ = (P[::-1] if revP else P), (Q[::-1] if revQ else Q)
                if swap:
                    PP, QQ = QQ, PP
                result.add(normalize_rows(cc, BB, TT+PP+QQ))
    return result


def regular_family(family):
    c, B, _ = family
    rows = [(c[:7], B[:7]), (c[7:], B[7:]),
            (tuple((c[j]-c[j+7]) % MOD for j in range(7)),
             tuple(tuple(x-y for x, y in zip(B[j], B[j+7])) for j in range(7)))]
    return all(oriented_cycles(cc, BB) for cc, BB in rows)


def maximal_families(families):
    """Remove repetitions and every component contained in a larger one."""
    maximal = []
    for d in (3, 2, 1, 0):
        layer = sorted(f for f in families if dimension(f) == d)
        previous = list(maximal)
        keep = [f for f in layer if not any(contains(g, f) for g in previous)]
        print(f'dimension {d}: {len(layer)} distinct -> {len(keep)} maximal', flush=True)
        maximal.extend(keep)
    return maximal


def quotient_maximal(maximal):
    left = set(maximal)
    classes = []
    while left:
        f = min(left, key=lambda f: (-dimension(f), f))
        orb = normalized_orbit(f)
        assert f in orb
        covered = orb & left
        classes.append((f, len(covered)))
        left -= covered
        print(f'class {len(classes)}: dimension {dimension(f)}, {len(covered)} normalized representatives', flush=True)
    return classes


def serialize(family):
    c, B, pivot = family
    return {'constant_24ths': c, 'basis': B, 'pivot': pivot}


if __name__ == '__main__':
    start = perf_counter()
    directory = Path(__file__).parent
    records = json.loads((directory/'cycle_components.json').read_text())['components']
    unique = {from_record(r) for r in records}
    print(f'{len(unique)} distinct normalized components', flush=True)
    maximal = maximal_families(unique)
    assert all(regular_family(f) for f in maximal)
    classes = quotient_maximal(maximal)
    result = {'distinct_components': len(unique),
              'maximal_counts': dict(Counter(dimension(f) for f in maximal)),
              'classes': [{'family': serialize(f), 'normalized_count': count} for f, count in classes]}
    (directory/'maximal_classes.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f'Elapsed: {perf_counter()-start:.2f}s', flush=True)
