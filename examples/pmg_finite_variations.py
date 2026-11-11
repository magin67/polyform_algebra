"""Finite variations of a unit triangle from PMG Variations, pmg-notes v5.

https://github.com/magin67/pmg-notes/tree/63d4a45da9f3a3cbe46f57f4e91c838f9a828c0c/pmg-library/en/Graph%20Transformations%20and%20Gluing
"""

from fractions import Fraction

from polyform_algebra.core import network, spanning_tree_weight
from polyform_algebra.variations import vary_resistances


def main():
    triangle = network(tuple("abc"), [("a", "b", 1), ("a", "c", 1), ("b", "c", 1)])
    updated = vary_resistances(triangle, {("a", "b"): Fraction(1, 3), ("a", "c"): Fraction(1, 3)})
    assert spanning_tree_weight(updated).value == Fraction(9, 5)
    print("Updated conductances:", [str(weight) for _, _, weight in updated.edges])


if __name__ == "__main__":
    main()
