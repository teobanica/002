"""Independent finite check: backtracking and direct row/column enumeration.

Does not import classify_ab.py or use its canonical invariant. It verifies all
25 candidates against the six literal representatives, visiting every one of
the 6*7! row/column permutations for each candidate. Phases cancel by dephasing.
"""
from itertools import permutations, combinations
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent


def search(p):
    found = []

    def extend(prefix, values, differences):
        j = len(prefix)
        if j == p:
            found.append(tuple(prefix))
            return
        for value in range(p):
            difference = (value-j) % p
            if value not in values and difference not in differences:
                extend(prefix+[value], values | {value}, differences | {difference})

    extend([0], {0}, {0})
    return found


def dephase(H, rows, columns, m):
    base, first, second = rows
    anchor = columns[0]
    return tuple((H[r][j]-H[base][j]-H[r][anchor]+H[base][anchor]) % m
                 for r in (first, second) for j in columns)


def is_five_plus_two(exponents):
    for five in combinations(range(7), 5):
        values = {(exponents[j]-exponents[five[0]]) % 20 for j in five}
        rest = [j for j in range(7) if j not in five]
        if values == {0, 4, 8, 12, 16} and \
                (exponents[rest[0]]-exponents[rest[1]]) % 20 == 10:
            return True
    return False


def bbb_sign(H):
    """Intrinsic triangle sign from the unique three antipodal pairs."""
    pairs = []
    for r, s in ((0, 1), (0, 2), (1, 2)):
        quotient = [(H[r][j]-H[s][j]) % 20 for j in range(7)]
        found = [{j, k} for j, k in combinations(range(7), 2)
                 if (quotient[j]-quotient[k]) % 20 == 10]
        assert len(found) == 1
        pairs.append(found[0])
    U, V, W = pairs
    assert len(U & V) == len(V & W) == 1
    k, = U & V
    j, = V & W
    return (H[1][k]-H[0][k]-H[1][j]+H[0][j]) % 20


def main():
    published = json.loads((ROOT/'classification.json').read_text())
    aa = search(7)
    bb = search(5)
    assert aa == [tuple(p) for p in published['AAA']['admissible']]
    assert bb == [tuple(p) for p in published['BBB']['admissible']]
    reps = {
        'A1': (7, [[0]*7, list(range(7)), [0, 2, 4, 6, 1, 3, 5]]),
        'A2': (7, [[0]*7, list(range(7)), [0, 3, 6, 2, 5, 1, 4]]),
        'A3': (7, [[0]*7, list(range(7)), [0, 2, 5, 1, 6, 4, 3]]),
        'A4': (7, [[0]*7, list(range(7)), [0, 2, 6, 5, 3, 1, 4]]),
        'B+': (20, [[0]*7, [0, 4, 8, 12, 16, 5, 15],
                    [15, 16, 12, 8, 4, 5, 0]]),
        'B-': (20, [[0]*7, [0, 4, 8, 12, 16, 15, 5],
                    [5, 16, 12, 8, 4, 15, 0]]),
    }
    targets = {7: {}, 20: {}}
    for label, (m, H) in reps.items():
        targets[m][dephase(H, (0, 1, 2), tuple(range(7)), m)] = label
    candidates = []
    for q in aa:
        expected = next(c['label'] for c in published['AAA']['classes']
                        if list(q) in c['members'])
        candidates.append((7, [[0]*7, list(range(7)), list(q)], expected,
                           {'type': 'AAA', 'pi': q}))
    for s in (1, -1):
        for q in bb:
            H = [[0]*7, [0, 4, 8, 12, 16, 5*s % 20, -5*s % 20],
                 [-5*s % 20]+[4*q[j] for j in range(1, 5)]+[5*s % 20, 0]]
            candidates.append((20, H, 'B+' if s == 1 else 'B-',
                               {'type': 'BBB', 'sign': s, 'pi': q}))
    certificates = []
    visited = 0
    column_orders = list(permutations(range(7)))
    for m, H, expected, description in candidates:
        for r, s in combinations(range(3), 2):
            quotient = [(H[r][j]-H[s][j]) % m for j in range(7)]
            assert (len(set(quotient)) == 7 if m == 7
                    else is_five_plus_two(quotient))
        matches = {}
        for rows in permutations(range(3)):
            # This also checks the triangle sign for all row reorderings.
            if m == 20:
                assert bbb_sign([H[r] for r in rows]) == (5 if expected == 'B+' else 15)
            for columns in column_orders:
                visited += 1
                value = dephase(H, rows, columns, m)
                label = targets[m].get(value)
                if label is not None:
                    matches.setdefault(label, {'rows': [r+1 for r in rows],
                                               'columns': [j+1 for j in columns]})
        assert set(matches) == {expected}, (description, expected, matches)
        certificate = matches[expected]
        assert dephase(H, tuple(r-1 for r in certificate['rows']),
                       tuple(j-1 for j in certificate['columns']), m) == \
            dephase(reps[expected][1], (0, 1, 2), tuple(range(7)), m)
        certificates.append({**description, 'class': expected, **certificate})
    # Verify each primary orbit certificate entry by entry, independently.
    for c in published['AAA']['classes']:
        H = [[0]*7, list(range(7)), c['representative']]
        for cert in c['certificates_from_representative']:
            transformed = dephase(H, tuple(r-1 for r in cert['rows']),
                                  tuple(j-1 for j in cert['columns']), 7)
            assert transformed == tuple(range(7))+tuple(cert['pi'])
    assert visited == 756000
    result = {'algorithm': 'Independent backtracking; exhaustive row/column permutations.',
              'visited_row_column_permutations': visited,
              'all_candidates_valid': True, 'all_primary_certificates_valid': True,
              'unique_target_class_for_each_candidate': True,
              'certificates_to_representatives': certificates}
    (ROOT/'independent_checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f'Independent check: {len(candidates)} candidates; {visited} row/column permutations.')
    print('Every candidate is valid and matches exactly its claimed class.')
    print('All covering certificates and all six BBB row-order sign checks passed.')


if __name__ == '__main__':
    main()
