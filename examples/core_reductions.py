"""Compare deletion, contraction and Kron reduction with exact arithmetic.

PMG v5, weighted path example and 'What should remain after elimination':
https://github.com/magin67/pmg-notes/blob/c909ddf8a590c69ed0a82fe48ca496e88fb707e0/pmg-library/en/Resistance%20Geometry%20and%20Duality/Spaces/Vertex%20Elimination%20and%20Induced%20Couplings.md
"""
from polyform_algebra.core import (
    network, dirichlet, extend, kron_reduce, contract, induced_subgraph,
    resistance, to_numeric,
)


def main():
    net = network(('a', 'x', 'y', 'b'),
                  (('a', 'x', 2), ('x', 'y', 1), ('y', 'b', 3)))
    d = dirichlet(net, ('a', 'b')).value
    print('Influence B:', d.B.values)
    print('Harmonic values:', extend(d, {'a': 0, 'b': 11}).value.values)
    graphs = {
        'Deletion': induced_subgraph(net, ('a', 'b')).value,
        'Contraction': contract(net, {'A': ('a', 'x'), 'B': ('y', 'b')}).value.network,
        'Kron reduction': kron_reduce(net, ('a', 'b')).value.network,
    }
    for name, graph in graphs.items():
        print(name, 'edges:', graph.edges, 'resistance:', resistance(graph).value.values[0][1])
    numeric = kron_reduce(to_numeric(net).value, ('a', 'b'))
    print('Numeric Kron edges:', numeric.value.network.edges)
    print('Grounded solve diagnostics:', numeric.diagnostics.numerical)


if __name__ == '__main__':
    main()
