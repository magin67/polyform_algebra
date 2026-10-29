"""Numerical examples from PMG introductory notes 6 and 7.

Sources are linked in the corresponding example modules. Expected constants
come from the notes; the independent matrix checks do not use polyform products.
"""
import numpy as np
import pytest
from examples.introductory.note06_single_edge import calculate as single_edge
from examples.introductory.note07_two_edges import calculate as two_edges


def test_note06_single_edge_values():
    """Tree weights 2+2+1 give tau=5; parallel paths give R12=2/5 and R13=3/5."""
    result=single_edge()
    expected={
        "tau_before":3, "R12_before":2/3, "R13_before":2/3, "R23_before":2/3,
        "coupling":1/3, "tau_after":5, "R12_after":2/5,
        "R13_after":3/5, "R23_after":3/5,
        "edge_probability_before":2/3, "edge_probability_after":4/5,
        "forest13_before":2, "forest13_after":3, "forest13_change":1,
        "tau_formula":5, "R13_formula":3/5, "forest13_change_formula":1,
    }
    for name,value in expected.items():
        assert result[name]==pytest.approx(value),name
    assert result["factor_update_matches_exp"]


def test_note07_two_edge_values_and_orientations():
    """Tree weights 4+2+2 give tau=8; two unit-resistance branches give R23=1/2."""
    result=two_edges()
    np.testing.assert_allclose(result["gram"],[[2/3,1/3],[1/3,2/3]])
    np.testing.assert_allclose(result["q"],[-1/3,1/3])
    expected={"tau_before":3,"tau_after_first":5,"R13_after_first":3/5,
              "tau_after":8,"R23_after":1/2,"tau_ratio_formula":8/3,
              "tau_formula":8,"R23_formula":1/2,"tau_sequential_formula":8,
              "mixed_coefficient":1,"mixed_coefficient_formula":1,
              "joint_probability_before":1/3}
    for name,value in expected.items():
        assert result[name]==pytest.approx(value),name
    assert result["factor_update_matches_exp"]
    assert result["reverse_order_matches"]


@pytest.mark.parametrize("calculate,weights", [(single_edge,(2,1,1)),(two_edges,(2,2,1))])
def test_variation_matches_independent_network_solution(calculate,weights):
    """Kirchhoff's cofactor and a grounded current solve independently check tau and R23."""
    c12,c13,c23=weights
    matrix=np.array([[c12+c13,-c12,-c13],[-c12,c12+c23,-c23],
                     [-c13,-c23,c13+c23]],dtype=float)
    result=calculate()
    assert result["tau_after"]==pytest.approx(np.linalg.det(matrix[:2,:2]))
    # Ground vertex 3; inject unit current at vertex 2 and remove it at vertex 3.
    potentials=np.linalg.solve(matrix[:2,:2],[0,1])
    assert result["R23_after"]==pytest.approx(potentials[1])
