"""Stage-two contexts and exact coefficients; PMG core specification v5, rev. 27.

https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Getting%20Started%20and%20Foundations/Foundations/Core%20specification.md
Sections: underlying space, parallel conductors, boundary operator.
"""
from fractions import Fraction
from polyform_algebra import Point, Simplex
from polyform_algebra.core import coordinates, network, reorder, to_numeric


def main():
    net = network(['a', 'b', 'isolated'], [('a', 'b', Fraction(1, 3)), ('b', 'a', Fraction(1, 6))])
    print('Aggregated conductance:', net.edges)
    print('Components:', net.components)
    vector = coordinates(net.context, [1, -1, 0]) / 3
    print('Exact coordinates:', vector.values)
    print('Reordered:', reorder(vector, ['isolated', 'b', 'a']).value.values)
    print('Numeric conversion:', to_numeric(net).diagnostics.status)
    a, b = Point.create_list('ab')
    print('Small exact boundary:', (Fraction(1, 10**12) * Simplex([a, b])).boundary())


if __name__ == '__main__':
    main()
