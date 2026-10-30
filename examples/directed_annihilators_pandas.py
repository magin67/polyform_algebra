"""Left and right Laplacian annihilators for three directed graphs.

Run: python -m examples.directed_annihilators_pandas
Optional dependency: python -m pip install pandas
PMG v5 revision d3b9cb30f8a990ebdc4c46b9cf4deb7f34d7a661:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/Balanced%20Graphs%20and%20Spanning-Tree%20Potentials.md
Sections 'Definition and Laplacian matrix' and 'Why positivity matters'.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/Balancing%20a%20Directed%20Graph%20by%20Spanning-Tree%20Potentials.md
Section 'A three-participant example'. Balancing is not performed here.
"""

import numpy as np

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis
from polyform_algebra.matrix.pandas import to_dataframe, to_series


def format_value(value):
    """Round the display only, avoiding confusing negative-zero output."""
    rounded = round(float(value), 6)
    return f"{0.0 if rounded == 0 else rounded:.6f}"


def main():
    try:
        import pandas as pd
    except ImportError:
        print("This example requires pandas: python -m pip install pandas")
        return

    cases = [
        ("Balanced directed cycle a -> b -> c -> a",
         [[1, -1, 0], [0, 1, -1], [-1, 0, 1]]),
        ("Unbalanced strongly connected graph (PMG example)",
         [[3, -1, -2], [-2, 3, -1], [-1, -1, 2]]),
        ("Two-dimensional kernels: a -> b and a -> c",
         [[2, -1, -1], [0, 0, 0], [0, 0, 0]]),
    ]
    for title, values in cases:
        axis = MatrixAxis(basis=Basis(["a", "b", "c"], space=Space(title, 3)))
        lap = Laplacian(values, rows=axis)
        result = lap.annihilators()
        print("\n" + title)
        print("Laplacian:")
        print(to_dataframe(lap).to_string(float_format=format_value))
        print("Kernel dimensions (left, right):", result.left_dimension, result.right_dimension)
        print("Left basis columns (U.T @ L = 0):")
        print(to_dataframe(result.left_basis).to_string(float_format=format_value))
        print("Right basis columns (L @ V = 0):")
        print(to_dataframe(result.right_basis).to_string(float_format=format_value))
        if result.left is not None:
            print("Unique annihilators normalized to sum n=3:")
            print(pd.concat([to_series(result.left), to_series(result.right)], axis=1)
                  .to_string(float_format=format_value))
        else:
            print("No unique normalized annihilator: left=None, right=None.")
            print("Left vectors: (0,s,t); right vectors: ((s+t)/2,s,t).")
        a = lap.to_numpy()
        print("Residual norms:", np.linalg.norm(result.left_basis.to_numpy().T @ a),
              np.linalg.norm(a @ result.right_basis.to_numpy()))

    print("\nThe unbalanced example has left=(15/16,12/16,21/16), right=(1,1,1).")
    print("Basis columns are Euclidean-coordinate orthonormal; a multidimensional basis is not canonical.")
    print("The graph weights have not been changed; balancing is a separate next step.")


if __name__ == "__main__":
    main()
