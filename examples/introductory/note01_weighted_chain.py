"""PMG notes 1-2: a three-vertex chain with conductances 1 and 2.

Sources (sections: chain Laplacian, potentials/power, three-vertex chain):
https://github.com/magin67/pmg-notes/blob/f63295490f9c4b35ef4ad418faa3cc75e5d1d0e7/pmg-library/en/Introduction/Laplacian%20-%20graph,%20electrical%20network%20and%20quadratic%20form.md
https://github.com/magin67/pmg-notes/blob/f63295490f9c4b35ef4ad418faa3cc75e5d1d0e7/pmg-library/en/Introduction/Effective%20Resistance%20and%20Graph%20Geometry.md
"""
import networkx as nx
import numpy as np

from polyform_algebra import Polyform, quForm


def calculate():
    """Build the network metric, then solve the note's 2 A current experiment."""
    graph = nx.Graph()
    # Weights are conductances, not edge resistances.
    graph.add_weighted_edges_from([(1, 2, 1), (2, 3, 2)])
    laplacian, points = Polyform.from_networkx(graph, return_points="dict")
    a, b, c = (points[label] for label in (1, 2, 3))
    f12, f23 = quForm([(a, b)]), quForm([(b, c)])
    metric = laplacian.exp()

    matrix, index = laplacian.to_matrix()
    # Respect the returned point-to-row mapping instead of assuming its order.
    order = [index[points[label]] for label in (1, 2, 3)]
    matrix = matrix[np.ix_(order, order)]
    # Ground vertex 3. This small NumPy solve is example code, not a network API.
    currents = np.array([2.0, 0.0, -2.0])
    potentials = np.append(np.linalg.solve(matrix[:2, :2], currents[:2]), 0.0)
    voltage12, voltage23 = potentials[0]-potentials[1], potentials[1]-potentials[2]
    return {
        "laplacian_matrix": matrix,
        "tau": metric.potential(quForm.one()),
        "R12": metric.norm(f12),
        "R23": metric.norm(f23),
        "R13": metric.norm(quForm([(a, c)])),
        # The only spanning tree uses both edges, with weight 1*2.
        "factorization_matches": metric.is_equivalent(
            (Polyform.one() + f12) * (Polyform.one() + 2*f23)
        ),
        "potentials": potentials,
        "edge_currents": np.array([voltage12, 2*voltage23]),
        "power": float(potentials @ matrix @ potentials),
        "edge_power": float(voltage12**2 + 2*voltage23**2),
        "source_power": float(potentials @ currents),
        "R13_from_current": float((potentials[0]-potentials[2])/2),
    }


def main():
    """Display the algebraic results and the independent electrical quantities."""
    print("PMG notes 1-2: weighted chain 1--2--3 (conductances 1, 2)")
    for name, value in calculate().items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
