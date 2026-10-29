"""Deletion examples: connected update, isolated vertices, and unnormalized potentials."""
import numpy as np
import pytest

from examples.introductory.note10_edge_deletion import calculate
from polyform_algebra import quForm


def test_cycle_edge_deletion_values():
    """Deleting 12 leaves the three-edge path 1-4-3-2, so tau=1 and R12=3."""
    result = calculate()["cycle"]
    assert result["tau_before"] == pytest.approx(4)
    assert result["removed_pair_R_before"] == pytest.approx(3/4)
    assert result["tau_after"] == pytest.approx(1)
    assert result["deletion_ratio"] == pytest.approx(1/4)
    assert result["removed_pair_R_after"] == pytest.approx(3)
    assert result["top_grade_after"] == 3


def test_bridge_deletion_preserves_vertex_set_and_rejects_global_norm():
    """Deleting conductance 2 isolates vertex 3; global tau=0 but the exponential survives."""
    result = calculate()["bridge"]
    updated, points = result["after"], result["points"]
    assert result["conductance"] == 2
    assert result["tau_before"] == pytest.approx(2)
    assert result["removed_pair_R_before"] == pytest.approx(1/2)
    assert result["tau_after"] == 0
    assert result["deletion_ratio"] == pytest.approx(0)
    assert set(updated.vertices) == set(points.values())
    assert result["top_grade_after"] == 1
    assert result["endpoint_potential_after"] == pytest.approx(1)
    assert result["removed_pair_R_after"] is None
    # Even a pair inside one component needs a separately constructed component metric.
    for pair in [(1, 2), (1, 3), (2, 3)]:
        with pytest.raises(ValueError, match="zero spanning coefficient"):
            updated.norm(quForm([(points[pair[0]], points[pair[1]])]))


@pytest.mark.parametrize("name", ["cycle", "bridge"])
def test_deletion_matches_rebuild_and_independent_laplacian(name):
    """Explicit remaining-edge matrices check tau and preserve the zero row of an isolated vertex."""
    result = calculate()[name]
    assert result["rebuild_matches"]
    assert result["restoration_matches"]
    n = 4 if name == "cycle" else 3
    edges = [(0, 3), (3, 2), (2, 1)] if name == "cycle" else [(0, 1)]
    matrix = np.zeros((n, n))
    for i, j in edges:
        matrix[i, i] += 1
        matrix[j, j] += 1
        matrix[i, j] -= 1
        matrix[j, i] -= 1
    actual, index = result["after"].extract_grade(1).to_matrix()
    order = [index[result["points"][label]] for label in range(1, n+1)]
    np.testing.assert_allclose(actual[np.ix_(order, order)], matrix)
    assert result["tau_after"] == pytest.approx(np.linalg.det(matrix[:-1, :-1]))
    if name == "cycle":
        # Ground vertex 4; inject at 1 and remove at 2.
        potential = np.linalg.solve(matrix[:-1, :-1], [1, -1, 0])
        assert result["removed_pair_R_after"] == pytest.approx(potential[0]-potential[1])
    else:
        p = result["points"]
        assert result["after"].is_equivalent(1*quForm.one()+quForm([(p[1], p[2])]))


@pytest.mark.parametrize("name", ["cycle", "bridge"])
def test_unnormalized_potential_update_survives_deletion(name):
    """Linearity of top-coefficient extraction holds even when deletion makes tau zero."""
    result = calculate()[name]
    before, after = result["before"], result["after"]
    p = result["removed_form"]
    points = result["points"]
    queries = [quForm.one(), p, quForm([(points[1], points[max(points)])])]
    for query in queries:
        expected = before.potential(query)-result["conductance"]*before.potential(p*query)
        assert after.potential(query) == pytest.approx(expected)
    assert result["tau_after"] == pytest.approx(result["tau_before"]*result["deletion_ratio"])
