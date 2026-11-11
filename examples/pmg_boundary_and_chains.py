"""Reproduce the two-path gluing and the square-ladder continuant in PMG v5.

https://github.com/magin67/pmg-notes/tree/4899950bb3620eb2294ad3d1a5872152d071015a/pmg-library/en/Graph%20Transformations%20and%20Gluing
"""

from polyform_algebra.boundary_data import forest_passport, pair_forest_passports
from polyform_algebra.core import network
from polyform_algebra.core_chains import cycle_chain_tree_count


def main():
    boundary = tuple("ab")
    first = network(tuple("abx"), [("a", "x", 1), ("x", "b", 1)])
    second = network(tuple("aby"), [("a", "y", 1), ("y", "b", 1)])
    direct = network(tuple("abxy"), first.edges + second.edges)
    combined = pair_forest_passports(
        boundary, forest_passport(first, boundary), forest_passport(second, boundary)
    )
    assert combined == forest_passport(direct, boundary)
    print("Two-path passport:", combined)
    count = cycle_chain_tree_count([4, 4, 4], [1, 1])
    assert count == 56
    print("Spanning trees of the three-square chain:", count)


if __name__ == "__main__":
    main()
