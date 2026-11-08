"""Boundary geometry of P5: the spanning-tree layer, selected vertices 0,2,4.

Source: PMG v5 4bea8e6e52c1299510dab05e5d96dbe6ea79d796,
https://github.com/magin67/pmg-notes/tree/v5/pmg-library/en/Graph%20Transformations%20and%20Gluing/Boundary%20polynomials%20and%20graph%20gluing
"""
from fractions import Fraction

from polyform_algebra.core import network, kron_reduce, resistance, spanning_tree_weight


def main():
    net = network(range(5), [(i, i+1, 1) for i in range(4)])
    reduced = kron_reduce(net, [0, 2, 4]).value
    # Each two-edge unit path becomes one edge with conductance 1/2.
    assert spanning_tree_weight(net).value == 1
    assert spanning_tree_weight(reduced.network).value == Fraction(1, 4)
    assert resistance(reduced.network).value.values[0][2] == 4
    print('P5 boundary (0,2,4): tau=1, reduced tau=1/4, R(0,4)=4')


if __name__ == '__main__':
    main()
