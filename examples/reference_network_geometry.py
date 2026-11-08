"""PMG stage-three exact and numeric geometry, matching public PythonExamples.

https://magin67.github.io/pmg-notes/ru/pythonexamples/
Pinned mathematical source: pmg-notes v5 c909ddf8a590c69ed0a82fe48ca496e88fb707e0.
"""
from polyform_algebra.core import (
    network, green, resistance, spanning_tree_weight, coordinates, gram,
    simplex_volume_squared, to_numeric, network_from_resistance,
)


def main():
    chain = network(['a','b','c'], [('a','b',1),('b','c',2)])
    print('Exact chain Green matrix:', green(chain).value.values)
    r = resistance(chain).value
    print('Exact chain resistances:', r.values)
    print('Reconstructed edges:', network_from_resistance(r.values,chain.context.labels).value.edges)
    triangle = network(['a','b','c'], [('a','b',1),('b','c',1),('a','c',1)])
    u = coordinates(triangle.context,[-1,1,0])
    v = coordinates(triangle.context,[-1,0,1])
    print('Triangle tree weight:', spanning_tree_weight(triangle).value)
    print('Triangle Gram matrix:', gram(triangle,[u,v]).value.values)
    print('Triangle squared area:', simplex_volume_squared(triangle,['a','b','c']).value)
    numeric = to_numeric(chain).value
    result = resistance(numeric,pairs=[('a','c')])
    print('Numeric endpoint resistance:',result.value.values)
    print('Numerical checks:',result.diagnostics.numerical)


if __name__ == '__main__':
    main()
