"""Exact calculations for the PMG Balanced Graphs series.

https://magin67.github.io/pmg-notes/en/directed-graphs-and-extensions/balanced-graphs/
"""

from polyform_algebra.directed import DirectedNetwork


def main():
    cycle = DirectedNetwork(tuple("abc"), (("a", "b", 1), ("b", "c", 1), ("c", "a", 1)))
    assert cycle.tree_potentials() == (1, 1, 1)
    assert cycle.forest_toponome() == {(): 1, (1,): 3, (2,): 3}
    print(cycle.one_way_resistances())


if __name__ == "__main__":
    main()
