"""Independent forest and matrix checks for the PMG note 10 examples."""
from itertools import combinations
from math import factorial

import networkx as nx
import numpy as np
import pytest

from examples.introductory.note10_path_cycle import calculate
from polyform_algebra import Polyform, quForm


@pytest.mark.parametrize("name,sums,tau", [
    ("P5", [1, 4, 6, 4, 1], 1),
    ("C4", [1, 4, 6, 4, 0], 4),
])
def test_path_cycle_grade_weights_and_top_form(name, sums, tau):
    """All path subsets are forests; only the complete four-cycle is excluded."""
    result = calculate()[name]
    assert result["forest_weight_sums"] == pytest.approx(sums)
    assert result["tau"] == pytest.approx(tau)
    points = list(result["points"].values())
    top = Polyform({quForm([tuple(points)]): tau}, vertices=points)
    assert result["grades"][len(points)-1].is_equivalent(top)
    assert result["factorization_matches"]
    # Four C4 trees combine into one stored top form with coefficient four.
    assert result["stored_term_counts"][len(points)-1] == 1


@pytest.mark.parametrize("name", ["P5", "C4"])
def test_grade_forms_against_independent_forest_enumeration(name):
    """Enumerate edge subsets and build component boundary forms without multiplying edges."""
    result = calculate()[name]
    graph = nx.path_graph(range(1, 6)) if name == "P5" else nx.cycle_graph(range(1, 5))
    points = result["points"]
    for k in range(5):
        terms = {}
        forest_count = 0
        for edges in combinations(graph.edges(), k):
            forest = nx.Graph()
            forest.add_nodes_from(graph.nodes())
            forest.add_edges_from(edges)
            if not nx.is_forest(forest):
                continue
            forest_count += 1
            # Singleton components contribute the scalar boundary 1.
            components = [tuple(points[i] for i in sorted(component))
                          for component in nx.connected_components(forest)
                          if len(component) > 1]
            form = quForm(components) if components else quForm.one()
            terms[form] = terms.get(form, 0) + 1
        expected = Polyform(terms, vertices=points.values())
        assert result["grades"][k].is_equivalent(expected)
        assert result["forest_weight_sums"][k] == forest_count
        # The factorial removes the k! orders of each distinct edge subset.
        power_grade = (1/factorial(k)) * (result["laplacian"]**k)
        assert result["grades"][k].is_equivalent(power_grade)


@pytest.mark.parametrize("name", ["P5", "C4"])
def test_all_pair_resistances_against_network_and_series_parallel_rules(name):
    """Kirchhoff's cofactor and grounded unit-current solves check every pair independently."""
    result = calculate()[name]
    n = 5 if name == "P5" else 4
    edges = [(i, i+1) for i in range(n-1)]
    if name == "C4":
        edges.append((n-1, 0))
    matrix = np.zeros((n, n))
    for i, j in edges:
        matrix[i, i] += 1
        matrix[j, j] += 1
        matrix[i, j] -= 1
        matrix[j, i] -= 1
    reduced = matrix[:-1, :-1]
    assert result["tau"] == pytest.approx(np.linalg.det(reduced))
    for (i, j), resistance in result["resistances"].items():
        current = np.zeros(n)
        current[i-1], current[j-1] = 1, -1
        potential = np.append(np.linalg.solve(reduced, current[:-1]), 0)
        assert resistance == pytest.approx(potential[i-1]-potential[j-1])
        distance = j-i
        expected = distance if name == "P5" else distance*(n-distance)/n
        assert resistance == pytest.approx(expected)
        assert result["potentials"][(i, j)] == pytest.approx(result["tau"]*expected)
