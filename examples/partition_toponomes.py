"""Forest and color-class examples from Partition Toponomes in PMG v5.

https://magin67.github.io/pmg-notes/en/polyform-algebra-and-combinatorics/partition-toponomes/
"""

from polyform_algebra.core import network
from polyform_algebra.toponomes import (
    coefficient_table, color_class_toponome, forest_toponome, partition_product,
)


def main():
    # From Polyforms to Partitions: the forest representation of the five-vertex path.
    path = network(tuple("abcde"), tuple((a, b, 1) for a, b in zip("abcd", "bcde")))
    forests = forest_toponome(path)
    assert coefficient_table(forests)[2, 1] == coefficient_table(forests)[2, 2] == 3
    print("P5:", forests)

    # Graph Colorings and Partition Products: the two-edge color-class model.
    disjoint = network(tuple("abcd"), (("a", "b", 1), ("c", "d", 1)))
    colors = color_class_toponome(disjoint)
    assert colors == partition_product({(1, 1): 1}, {(1, 1): 1})
    print("Two disjoint edges:", colors)


if __name__ == "__main__":
    main()
