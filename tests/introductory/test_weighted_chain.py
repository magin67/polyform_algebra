"""Check the chain example from PMG notes 1-2 (sources in its module)."""
import numpy as np
import pytest

from examples.introductory.note01_weighted_chain import calculate


def test_chain_laplacian_and_metric():
    """Series resistances add: 1 + 1/2 = 3/2; the unique tree has weight 1*2."""
    result = calculate()
    np.testing.assert_allclose(
        result["laplacian_matrix"], [[1, -1, 0], [-1, 3, -2], [0, -2, 2]]
    )
    for name, expected in {"tau": 2, "R12": 1, "R23": 0.5, "R13": 1.5}.items():
        assert result[name] == pytest.approx(expected), name
    assert result["factorization_matches"]


def test_chain_against_independent_current_solutions():
    """A cofactor and three unit-current solves check the polyform evaluations."""
    # Construct the reference matrix explicitly, without library conversion.
    matrix = np.array([[1, -1, 0], [-1, 3, -2], [0, -2, 2]], dtype=float)
    result = calculate()
    assert result["tau"] == pytest.approx(np.linalg.det(matrix[:2, :2]))
    for i, j, name in [(0, 1, "R12"), (1, 2, "R23"), (0, 2, "R13")]:
        current = np.zeros(3)
        current[i], current[j] = 1, -1
        potential = np.append(np.linalg.solve(matrix[:2, :2], current[:2]), 0)
        assert result[name] == pytest.approx(potential[i]-potential[j])


def test_chain_current_and_power_balance():
    """With 2 A from 1 to 3, potentials are (3,1,0), edge currents 2 A, power 6 W."""
    result = calculate()
    np.testing.assert_allclose(result["potentials"], [3, 1, 0])
    np.testing.assert_allclose(result["edge_currents"], [2, 2])
    np.testing.assert_allclose(
        result["laplacian_matrix"] @ result["potentials"], [2, 0, -2]
    )
    for name in ("power", "edge_power", "source_power"):
        assert result[name] == pytest.approx(6), name
    assert result["power"] == pytest.approx(2**2 * result["R13"])
    assert result["R13_from_current"] == pytest.approx(result["R13"])
