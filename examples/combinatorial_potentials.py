"""Examples from the Combinatorial Potentials series in PMG.

https://magin67.github.io/pmg-notes/en/polyform-algebra-and-combinatorics/combinatorial-potentials/
"""

from math import factorial, perm

from polyform_algebra.combinatorics import (
    first_kind_potential,
    path_cycle_potential,
    second_kind_potential,
)


def main():
    # Partitions into two pairs, then permutations with two triples and a pair.
    assert second_kind_potential((2, 2)) == 3
    assert first_kind_potential((3, 3, 2)) == 1120

    # A directed path of two arcs on a, b, c: states and arc-addition histories.
    potential = path_cycle_potential(path_lengths=(2,))
    histories = potential * perm(len("abc"), 3)
    assert potential == 2 and histories // factorial(2) == 6 and histories == 12
    print("s₂([2,2]) = 3 | s₁([3,3,2]) = 1120 | u = 2, N₃ = 6, H₃ = 12")


if __name__ == "__main__":
    main()
