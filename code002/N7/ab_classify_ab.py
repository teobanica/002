"""Exact AAA/BBB classification after the reductions proved in the PDF.

Python 3 standard library only. No floating point, root-order search, or
identification by complex conjugation is used. Run: python3 run_ab.py
"""
from itertools import permutations
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
AAA_REPS = (
    (0, 2, 4, 6, 1, 3, 5),
    (0, 3, 6, 2, 5, 1, 4),
    (0, 2, 5, 1, 6, 4, 3),
    (0, 2, 6, 5, 3, 1, 4),
)


def admissible(p):
    """All permutations pi fixing zero with pi-id also a permutation."""
    answer = []
    for tail in permutations(range(1, p)):
        q = (0,) + tail
        if len({(q[j] - j) % p for j in range(p)}) == p:
            answer.append(q)
    return answer


def canonical(H, m):
    """Complete monomial invariant for a matrix of m-th roots.

    Choose the new first row, the order of the other two rows, and a column
    anchor; dephase, then sort the column pairs. Arbitrary unit phases cancel
    in dephasing. Thus six row orders times seven anchors exhaust all choices.
    """
    candidates = []
    for b, r, s in permutations(range(3)):
        columns = [((H[r][j] - H[b][j]) % m,
                    (H[s][j] - H[b][j]) % m) for j in range(7)]
        for x0, y0 in columns:
            candidates.append(tuple(sorted(((x-x0) % m, (y-y0) % m)
                                           for x, y in columns)))
    return min(candidates)


def aaa_matrix(pi):
    return [[0]*7, list(range(7)), list(pi)]


def bbb_matrix(sign, pi):
    return [[0]*7, [0, 4, 8, 12, 16, (5*sign) % 20, (-5*sign) % 20],
            [(-5*sign) % 20] + [4*pi[j] for j in range(1, 5)]
            + [(5*sign) % 20, 0]]


def aaa_orbit(pi):
    """All normalized permutations in one monomial orbit, with certificates."""
    H = aaa_matrix(pi)
    result = {}
    for b, r, s in permutations(range(3)):
        f = [(H[r][j]-H[b][j]) % 7 for j in range(7)]
        g = [(H[s][j]-H[b][j]) % 7 for j in range(7)]
        for anchor in range(7):
            order = sorted(range(7), key=lambda j: (f[j]-f[anchor]) % 7)
            q = tuple((g[j]-g[anchor]) % 7 for j in order)
            result.setdefault(q, {'rows': [b+1, r+1, s+1],
                                  'columns': [j+1 for j in order]})
    return result


def main():
    a = admissible(7)
    b = admissible(5)
    assert len(a) == 19 and len(b) == 3
    assert b == [(0, 2, 4, 1, 3), (0, 3, 1, 4, 2), (0, 4, 3, 2, 1)]
    a_keys = [canonical(aaa_matrix(q), 7) for q in AAA_REPS]
    assert len(set(a_keys)) == 4
    a_classes = [[q for q in a if canonical(aaa_matrix(q), 7) == key]
                 for key in a_keys]
    assert list(map(len, a_classes)) == [3, 2, 7, 7]
    orbits = [aaa_orbit(q) for q in AAA_REPS]
    for orbit, members in zip(orbits, a_classes):
        assert set(orbit) == set(members)
    assert set().union(*(set(o) for o in orbits)) == set(a)
    b_reps = [bbb_matrix(s, b[-1]) for s in (1, -1)]
    b_keys = [canonical(H, 20) for H in b_reps]
    assert b_keys[0] != b_keys[1]
    b_candidates = []
    for sign in (1, -1):
        for q in b:
            H = bbb_matrix(sign, q)
            target = b_keys.index(canonical(H, 20))
            assert target == (0 if sign == 1 else 1)
            b_candidates.append({'sign': sign, 'pi': q,
                                 'class': 'B+' if target == 0 else 'B-'})
    for i, q in enumerate(AAA_REPS):
        conjugate = [[(-x) % 7 for x in row] for row in aaa_matrix(q)]
        assert a_keys.index(canonical(conjugate, 7)) == [0, 1, 3, 2][i]
    result = {
        'conventions': 'Ordinary monomial equivalence; conjugation is not included.',
        'AAA': {'tested_permutations': 720, 'admissible': a,
                'class_sizes': list(map(len, a_classes)),
                'classes': [{'label': f'A{i+1}', 'representative': q,
                             'members': a_classes[i],
                             'certificates_from_representative': [
                                 {'pi': p, **orbits[i][p]} for p in sorted(orbits[i])]}
                            for i, q in enumerate(AAA_REPS)]},
        'BBB': {'tested_permutations': 24, 'admissible': b,
                'candidates': b_candidates,
                'representatives': [{'label': label, 'modulus': 20, 'exponents': H}
                                    for label, H in zip(('B+', 'B-'), b_reps)]},
    }
    (ROOT/'classification.json').write_text(json.dumps(result, indent=2)+'\n')
    print('AAA: 720 permutations; 19 admissible; 4 classes of sizes 3, 2, 7, 7.')
    print('BBB: 24 permutations; 3 admissible; 6 signed candidates; 2 classes.')
    print('Conjugation fixes A1,A2 and interchanges A3,A4 and B+,B-.')


if __name__ == '__main__':
    main()
