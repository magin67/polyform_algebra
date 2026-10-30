"""Compare circumspheres, spanning-tree and eigenvector centralities of P5 and K1,4.

Run: python -m examples.circumsphere_centrality_pandas
Optional dependency: python -m pip install pandas
Source: PMG v5 revision d3b9cb30f8a990ebdc4c46b9cf4deb7f34d7a661,
'Centroid, Spanning-Tree, and Eigenvector Centralities', equations (4)-(5),
and 'Eigenvector centrality and a common normalization',
'The complete graph, a path, and a star'.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Metric%20parameters/Centroid,%20Spanning-Tree,%20and%20Eigenvector%20Centralities.md
"""

import networkx as nx

from polyform_algebra.matrix.networkx import from_networkx
from polyform_algebra.matrix.pandas import to_series


def format_value(value):
    """Round only the display, suppressing negative zero from floating-point noise."""
    rounded = round(float(value), 6)
    return f"{0.0 if rounded == 0 else rounded:.6f}"


def main():
    try:
        import pandas as pd
    except ImportError:
        print("This example requires pandas: python -m pip install pandas")
        return

    for name, graph in [("Path P5: 0-1-2-3-4", nx.path_graph(5)),
                        ("Star K1,4: center 0, leaves 1,2,3,4", nx.star_graph(4))]:
        lap = from_networkx(graph)
        sphere = lap.circumsphere()
        centrality = lap.spanning_tree_centrality()
        eigen = lap.eigenvector_centrality()
        frame = pd.concat([
            to_series(sphere.center_coordinates),
            to_series(centrality.expected_degrees),
            to_series(centrality.scores),
            to_series(eigen.scores),
        ], axis=1)
        print(name)
        print(f"Resistance radius (squared linear radius): {sphere.resistance_radius:.6f}")
        print(f"Largest adjacency eigenvalue: {eigen.eigenvalue:.6f}")
        print(frame.to_string(float_format=format_value))
        print()

    # A unit tree has rho=sum(edge resistances)/4=1 in both examples.
    # Its unique spanning tree has the original unweighted vertex degrees.
    print("Both radii are 1; their center coordinates and centralities differ.")
    print("Path center: (vertex 0 + vertex 4)/2.")
    print("Star center: -vertex 0 + (vertex 1 + vertex 2 + vertex 3 + vertex 4)/2.")
    print("The negative center coordinate is valid: the center lies outside the simplex.")
    print("Expected tree degrees sum to 8, normalized centralities to 1.")
    print("Each center coordinate equals 1 - expected_tree_degree/2.")
    # P5's Perron vector is proportional to (1,sqrt(3),2,sqrt(3),1).
    # The star's vector is proportional to (2,1,1,1,1), with eigenvalue 2.
    print("On the path, eigenvector centrality singles out the middle vertex;")
    print("spanning-tree centrality assigns 1/4 to each of the three internal vertices.")
    print("At the star center, eigenvector centrality is 1/3, spanning-tree centrality is 1/2.")

    # Explicit component mode never introduces a finite metric between parts.
    disconnected = nx.Graph()
    disconnected.add_edge("a", "b", weight=2)
    disconnected.add_node("isolated")
    lap = from_networkx(disconnected)
    print("\nDisconnected example: one weighted edge and an isolated vertex")
    for sphere, centrality, eigen in zip(lap.circumsphere(componentwise=True),
                                        lap.spanning_tree_centrality(componentwise=True),
                                        lap.eigenvector_centrality(componentwise=True)):
        print("Component:", sphere.center_coordinates.labels)
        print("Resistance radius:", sphere.resistance_radius)
        if centrality.scores is None:
            print("Expected degree is 0; normalized centrality is undefined (None).")
        else:
            print(to_series(centrality.scores).to_frame().to_string())
        print(to_series(eigen.scores).to_frame().to_string())


if __name__ == "__main__":
    main()
