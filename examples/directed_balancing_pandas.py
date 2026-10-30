"""Balance a directed graph by its left Laplacian annihilator; display pandas tables.

Run: python -m examples.directed_balancing_pandas
Optional dependency: python -m pip install pandas
PMG v5 revision e28771c349cb54d84040702e7dc0cd345b57fd0f:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/Balancing%20a%20Directed%20Graph%20by%20Spanning-Tree%20Potentials.md
Equations (1)-(4), 'A three-participant example', and normalization (2).
"""

import numpy as np

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis
from polyform_algebra.matrix.pandas import to_dataframe, to_series


def format_value(value):
    """Suppress negative zero in display only; computations retain their precision."""
    rounded = round(float(value), 6)
    return f"{0.0 if rounded == 0 else rounded:.6f}"


def main():
    try:
        import pandas as pd
    except ImportError:
        print("This example requires pandas: python -m pip install pandas")
        return

    axis = MatrixAxis(basis=Basis(["a", "b", "c"], space=Space("participants", 3)))
    original = Laplacian([[3, -1, -2], [-2, 3, -1], [-1, -1, 2]], rows=axis)
    result = original.balance()
    balanced = result.laplacian
    print("Source multipliers (mean one):")
    print(to_series(result.multipliers).to_frame().to_string())
    print("Expected: (15/16, 12/16, 21/16).")

    for label, lap in [("Original", original), ("Balanced", balanced)]:
        adjacency = to_dataframe(lap.to_adjacency())
        print(f"\n{label} adjacency (rows are sources):")
        print(adjacency.to_string())
        print(f"{label} Laplacian:")
        print(to_dataframe(lap).to_string())
        totals = pd.DataFrame({"outgoing": adjacency.sum(axis=1),
                               "incoming": adjacency.sum(axis=0)})
        totals["outgoing_minus_incoming"] = totals["outgoing"] - totals["incoming"]
        print("Balance at each vertex:")
        print(totals.to_string(float_format=format_value))

    a = original.to_adjacency().to_numpy()
    b = balanced.to_adjacency().to_numpy()
    np.testing.assert_allclose(b.sum(axis=1), b.sum(axis=0), atol=1e-12)
    np.testing.assert_array_equal(a > 0, b > 0)
    np.testing.assert_allclose(a / a.sum(axis=1)[:, None], b / b.sum(axis=1)[:, None])
    np.testing.assert_allclose(b / result.multipliers.to_numpy()[:, None], a)
    assert balanced.vector_type() == "vector"
    print("\nChecks passed: balance, arc preservation, outgoing proportions, recovery.")
    print("Expected balanced incoming/outgoing weights: (45/16, 36/16, 42/16).")
    print("Mean-one multipliers do not preserve total flow:", a.sum(), "->", b.sum())
    print("The source graph is unchanged. The result is balanced, not necessarily symmetric.")


if __name__ == "__main__":
    main()
