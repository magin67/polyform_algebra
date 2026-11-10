"""Compute the degree potential and compare the two graph circumspheres.

PMG v5, Sphere Centers and Graph Heterogeneity, exact examples:
https://github.com/magin67/pmg-notes/blob/de21abc72638400e642cb5559af19d833a16f2e6/pmg-library/en/Resistance%20Geometry%20and%20Duality/Metric%20parameters/Sphere%20Centers%20and%20Graph%20Heterogeneity.md
"""

from polyform_algebra.core import circumsphere, laplacian_circumsphere, network


def main():
    net = network(('a', 'b', 'c'), (('a', 'b', 7), ('a', 'c', 1), ('b', 'c', 1)))
    laplacian = laplacian_circumsphere(net).value
    resistance = circumsphere(net).value

    print('Vertex order:', ', '.join(net.context.labels))
    print('Degree potential:', ', '.join(map(str, laplacian.degree_potential.values)))
    print('Laplacian center:', ', '.join(map(str, laplacian.center_coordinates.values)))
    print('Laplacian squared radius:', laplacian.squared_radius)
    print('Resistance center:', ', '.join(map(str, resistance.center_coordinates.values)))
    print('Resistance radius:', resistance.resistance_radius)


if __name__ == '__main__':
    main()
