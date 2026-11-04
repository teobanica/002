"""Support-pattern reduction for regular 3x7 partial Hadamard matrices.

A row with cycle type 3+2+2 determines a partition of seven columns into
one triple and two unordered pairs.  Fixing the first row's support partition,
we enumerate all supports available to a second row modulo the stabilizer of
the first partition.  This is the finite combinatorial front end to the exact
phase equations; it does not yet solve those equations.
"""
from itertools import combinations, permutations

COLUMNS = frozenset(range(7))
BASE = ((0, 1, 2), (3, 4), (5, 6))


def support_partitions():
    """All 105 partitions T|A|B with |T|=3 and unordered |A|=|B|=2."""
    result = []
    for triple in combinations(range(7), 3):
        rest = sorted(COLUMNS.difference(triple))
        for pair_a in combinations(rest, 2):
            pair_b = tuple(x for x in rest if x not in pair_a)
            if pair_a < pair_b:
                result.append((triple, pair_a, pair_b))
    return result


def base_stabilizer():
    """Permutations preserving BASE as a set of blocks."""
    result = set()
    triple, pair_a, pair_b = BASE
    for image_t in permutations(triple):
        for image_a in permutations(pair_a):
            for image_b in permutations(pair_b):
                for swap_pairs in (False, True):
                    target_a, target_b = ((image_a, image_b) if not swap_pairs
                                          else (image_b, image_a))
                    p = [None] * 7
                    for old, new in zip(triple, image_t):
                        p[old] = new
                    for old, new in zip(pair_a, target_a):
                        p[old] = new
                    for old, new in zip(pair_b, target_b):
                        p[old] = new
                    result.add(tuple(p))
    return sorted(result)


def act(partition, permutation):
    blocks = [tuple(sorted(permutation[i] for i in block))
              for block in partition]
    triple = next(block for block in blocks if len(block) == 3)
    pairs = sorted(block for block in blocks if len(block) == 2)
    return (triple, *pairs)


def support_orbits():
    """Return canonical orbit representatives and their sizes."""
    remaining = set(support_partitions())
    group = base_stabilizer()
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {act(representative, p) for p in group} & remaining
        orbits.append((representative, len(orbit)))
        remaining.difference_update(orbit)
    return orbits


def main():
    partitions = support_partitions()
    group = base_stabilizer()
    orbits = support_orbits()
    assert len(partitions) == 105
    assert sum(size for _, size in orbits) == len(partitions)
    print(f"support partitions: {len(partitions)}")
    print(f"stabilizer size: {len(group)}")
    print(f"support orbits: {len(orbits)}")
    for index, (representative, size) in enumerate(orbits, start=1):
        print(index, size, representative)


if __name__ == "__main__":
    main()
