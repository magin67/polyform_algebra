"""PMG weighted-chain transformations and disconnected resistance.

Source: PMG Introduction, 'Laplacian, Green Matrix, and Effective Resistance
Matrix', section 'Example: a three-vertex path'. Revision:
c27d4f781640c3639df8c0d6f1a9ca11795b8139, branch v5.
https://github.com/magin67/pmg-notes/tree/v5/pmg-library/en/Introduction
"""

import networkx as nx
from polyform_algebra.matrix.networkx import from_networkx


def main():
    """Conductances 1,2 yield series resistances 1,1/2,3/2."""
    graph = nx.Graph()
    graph.add_edge("a", "b", weight=1)
    graph.add_edge("b", "c", weight=2)
    laplacian = from_networkx(graph, vertices=["a", "b", "c"])
    print("Laplacian:", laplacian.vector_type())
    print(laplacian.to_numpy())
    greenian = laplacian.to_greenian()
    resistance = greenian.to_resistance()
    print("Greenian:")
    print(greenian.to_numpy())
    print("Resistance:")
    print(resistance.to_numpy())
    print("Recovered Laplacian:")
    print(resistance.to_laplacian().to_numpy())
    graph.add_node("isolated")
    disconnected = from_networkx(graph).to_resistance()
    print("Components:", disconnected.components)
    print(disconnected.to_numpy())


if __name__ == "__main__":
    main()
