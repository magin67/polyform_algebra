"""Compare the resistance circumspheres of the four-vertex path and star.

PMG v5, Sphere Centers and Graph Heterogeneity, exact examples:
https://github.com/magin67/pmg-notes/blob/de21abc72638400e642cb5559af19d833a16f2e6/pmg-library/en/Resistance%20Geometry%20and%20Duality/Metric%20parameters/Sphere%20Centers%20and%20Graph%20Heterogeneity.md
"""

from polyform_algebra.core import circumsphere, network


def main():
    vertices = ('a', 'b', 'c', 'd')
    graphs = (
        ('Path', (('a', 'b', 1), ('b', 'c', 1), ('c', 'd', 1))),
        ('Star', (('a', 'b', 1), ('a', 'c', 1), ('a', 'd', 1))),
    )
    for name, edges in graphs:
        sphere = circumsphere(network(vertices, edges)).value
        center = ', '.join(
            f'{vertex}={value}'
            for vertex, value in zip(vertices, sphere.center_coordinates.values)
        )
        print(f'{name}: radius = {sphere.resistance_radius}, center = ({center})')


if __name__ == '__main__':
    main()
