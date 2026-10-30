"""Compare induced deletion, contraction and elimination on the same path.

PMG v5, read 2026-10-30:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Vertex%20Elimination%20and%20Induced%20Couplings.md
Revision 5, 'What should remain after elimination': path weights 2, 1, 3.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Vertex%20Contraction%20and%20Graph%20Gluing.md
Revision 6 (Russian revision 4), vertex identification, equation (1).
Run: python -m examples.graph_operations_pandas
"""

import numpy as np

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis
from polyform_algebra.matrix.pandas import to_dataframe


def main():
    try:
        import pandas as pd
    except ImportError:
        print("This example requires pandas: python -m pip install pandas")
        return

    axis = MatrixAxis(basis=Basis(["a", "x", "y", "b"], space=Space("path", 4)))
    original = Laplacian([[2, -2, 0, 0], [-2, 3, -1, 0],
                          [0, -1, 4, -3], [0, 0, -3, 3]], rows=axis)
    print("Original path: a --(2)-- x --(1)-- y --(3)-- b")
    print("\nOriginal adjacency:")
    print(to_dataframe(original.to_adjacency()))
    print("\nOriginal Laplacian:")
    print(to_dataframe(original))

    induced = original.induced_subgraph(keep=["a", "b"])
    contraction = original.contract({"a": ["a", "x"], "b": ["y", "b"]})
    elimination = original.eliminate(keep=["a", "b"])
    graphs = {"Induced subgraph": induced,
              "Contraction": contraction.laplacian,
              "Elimination": elimination.laplacian}
    print("\nContraction membership: new a={a,x}, new b={y,b}")
    print(to_dataframe(contraction.membership_matrix))
    print("\nElimination: bi-coordinates of projections of x,y onto a,b")
    print(to_dataframe(elimination.dirichlet.influence_matrix).round(6))

    rows = []
    for name, lap in graphs.items():
        print(f"\n{name}: Laplacian")
        print(to_dataframe(lap).round(6))
        print("Adjacency:")
        print(to_dataframe(lap.to_adjacency()).round(6))
        resistance = to_dataframe(lap.to_resistance())
        print("Resistance:")
        print(resistance.round(6))
        rows.append({"operation": name,
                     "coupling_a_b": lap.to_adjacency().to_numpy()[0, 1],
                     "resistance_a_b": resistance.loc["a", "b"]})

    print("\nComparison (contraction labels denote groups, not the original points):")
    print(pd.DataFrame(rows).set_index("operation").round(6))
    # Series resistances add: 1/2+1+1/3=11/6. Only elimination preserves it.
    source_resistance = original.to_resistance().to_numpy()[0, 3]
    np.testing.assert_allclose(source_resistance, 11/6)
    np.testing.assert_allclose([row["coupling_a_b"] for row in rows], [0, 1, 6/11])
    assert np.isposinf(rows[0]["resistance_a_b"])
    np.testing.assert_allclose(rows[1]["resistance_a_b"], 1)
    np.testing.assert_allclose(rows[2]["resistance_a_b"], source_resistance)
    s = contraction.membership_matrix.to_numpy()
    np.testing.assert_allclose(contraction.laplacian.to_numpy(), s.T @ original.to_numpy() @ s)
    print("\nChecks passed. Original R(a,b)=11/6; only elimination retains that metric.")
    print("Deletion removes incident edges. Contraction identifies vertices. Elimination induces replacement couplings.")


if __name__ == "__main__":
    main()
