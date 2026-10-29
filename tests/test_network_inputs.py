"""Network input contracts; signed coefficients remain valid in the algebra.

Source: PMG core specification, scope and graph/Laplacian representation:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/ru/Основания/Спецификация%20ядра%20PMG.md
"""
import networkx as nx
import numpy as np
import pytest

from polyform_algebra import Point, Vector, Polyform, quForm


@pytest.mark.parametrize("matrix", [1, [1, -1], [[[0]]], [[1, 2]], [[1], [2, 3]]])
def test_laplacian_requires_square_two_dimensional_input(matrix):
    """Malformed shapes fail explicitly instead of indexing a nonexistent dimension."""
    with pytest.raises(ValueError, match="Matrix must"):
        Polyform.from_laplacian(matrix, Point.create_list("ab"))


@pytest.mark.parametrize("matrix", [[[np.nan]], [[np.inf]], [[-np.inf]]])
def test_laplacian_rejects_nonfinite_entries(matrix):
    """NaN and infinities cannot define finite network conductances."""
    with pytest.raises(ValueError, match="finite"):
        Polyform.from_laplacian(matrix, [Point("a")])


@pytest.mark.parametrize("matrix", [[["0"]], [[0j]], [[True]], [[None]]])
def test_laplacian_rejects_nonreal_numeric_input(matrix):
    """Strings, booleans, complex values and missing values are not silently coerced."""
    with pytest.raises(TypeError, match="real numbers"):
        Polyform.from_laplacian(matrix, [Point("a")])


@pytest.mark.parametrize("matrix,message", [
    ([[7,-1],[-1,7]], "row sums"),
    ([[1,-1],[-2,2]], "symmetric"),
    ([[-1,1],[1,-1]], "non-negative"),
    ([[0,-1],[-1,0]], "row sums"),
    ([[-1,0],[0,0]], "non-negative"),
])
def test_laplacian_rejects_invalid_network_matrices(matrix, message):
    """Do not discard a wrong diagonal or accept asymmetry and negative conductances."""
    with pytest.raises(ValueError, match=message):
        Polyform.from_laplacian(matrix, Point.create_list("ab"))


@pytest.mark.parametrize("tolerance", [-1, np.nan, np.inf, True, "small"])
def test_laplacian_rejects_invalid_tolerance(tolerance):
    """The absolute error bound must be a finite non-negative real number."""
    with pytest.raises(ValueError, match="Tolerance"):
        Polyform.from_laplacian([[0]], [Point("a")], tolerance=tolerance)


def test_laplacian_symmetry_check_has_no_relative_tolerance():
    """Large weights must not hide a discrepancy of one via NumPy's default rtol."""
    matrix = [[1e8,-1e8],[-1e8+1,1e8-1]]
    with pytest.raises(ValueError, match="symmetric"):
        Polyform.from_laplacian(matrix, Point.create_list("ab"))


def test_laplacian_averages_small_asymmetry_without_mutating_input():
    """Errors inside the absolute bound are symmetrized, retaining the mean weight."""
    delta = 1e-8
    matrix = np.array([[1,-1],[-1-delta,1+delta]])
    original = matrix.copy()
    p = Polyform.from_laplacian(matrix, Point.create_list("ab"), tolerance=2*delta)
    expected = (1+delta/2)*np.array([[1,-1],[-1,1]])
    np.testing.assert_allclose(p.to_matrix()[0], expected, atol=1e-12, rtol=0)
    np.testing.assert_array_equal(matrix, original)


@pytest.mark.parametrize("labels,error", [
    ([Point("a")], ValueError),
    ([Point("a"),Point("a")], ValueError),
    (["a","b"], TypeError),
    ([Point("a"),Vector("b")], TypeError),
])
def test_laplacian_requires_distinct_points_in_matrix_order(labels, error):
    """Mismatched, repeated or non-point labels must not merge or misidentify vertices."""
    with pytest.raises(error):
        Polyform.from_laplacian([[1,-1],[-1,1]], labels)


@pytest.mark.parametrize("size", [0,1,3])
def test_empty_networks_keep_all_vertices(size):
    """Zero Laplacians, including the empty matrix, retain their vertex dimensions."""
    labels = Point.create_list([str(i) for i in range(size)])
    for p in [Polyform.from_laplacian(np.zeros((size,size)), labels),
              Polyform.from_networkx(nx.empty_graph(size))]:
        assert p.is_zero()
        assert len(p.vertices) == size
        assert p.to_matrix()[0].shape == (size,size)
        assert p.exp().is_one()


def test_laplacian_round_trip_preserves_order_and_disconnected_components():
    """Two weighted components and one isolated vertex retain the supplied order."""
    labels = Point.create_list("edcba")
    matrix = np.array([[2,-2,0,0,0],[-2,2,0,0,0],
                       [0,0,3,-3,0],[0,0,-3,3,0],[0,0,0,0,0]])
    p = Polyform.from_laplacian(matrix, labels)
    actual,index = p.to_matrix()
    assert list(index) == labels
    np.testing.assert_allclose(actual, matrix)


@pytest.mark.parametrize("graph", [nx.DiGraph([(0,1)]), nx.MultiDiGraph([(0,1)])])
def test_directed_graphs_are_rejected(graph):
    """Import must not silently erase edge directions."""
    with pytest.raises(ValueError, match="Directed"):
        Polyform.from_networkx(graph)


@pytest.mark.parametrize("cls", [nx.Graph,nx.MultiGraph])
@pytest.mark.parametrize("weight", [0,1])
def test_self_loops_are_rejected_even_with_zero_weight(cls, weight):
    """A loop is outside the network model and must not silently disappear."""
    graph=cls()
    graph.add_edge(0,0,weight=weight)
    with pytest.raises(ValueError, match="Self-loops"):
        Polyform.from_networkx(graph)


@pytest.mark.parametrize("weight,error", [(-1,ValueError),(np.nan,ValueError),
    (np.inf,ValueError),(-np.inf,ValueError),("2",TypeError),
    (1j,TypeError),(True,TypeError),(None,TypeError)])
def test_graph_weights_must_be_finite_nonnegative_reals(weight,error):
    """Only physical conductances are accepted by the network importer."""
    graph=nx.Graph()
    graph.add_edge(0,1,weight=weight)
    with pytest.raises(error):
        Polyform.from_networkx(graph)


def test_parallel_conductances_are_summed_and_zero_edges_keep_vertices():
    """Parallel weights 2+3+1 give conductance 6; zero edges add no connection."""
    graph=nx.MultiGraph()
    graph.add_edge(0,1,weight=2)
    graph.add_edge(0,1,weight=np.float64(3))
    graph.add_edge(0,1)  # Missing weight means unit conductance.
    graph.add_edge(1,2,weight=0)
    p,points=Polyform.from_networkx(graph,return_points="list")
    np.testing.assert_array_equal(p.to_matrix()[0],[[6,-6,0],[-6,6,0],[0,0,0]])
    assert len(points)==3
    assert graph.number_of_edges()==4


def test_negative_parallel_edge_cannot_be_hidden_by_positive_sum():
    """Validate individual conductances before summing parallel connections."""
    graph=nx.MultiGraph()
    graph.add_edge(0,1,weight=2)
    graph.add_edge(0,1,weight=-1)
    with pytest.raises(ValueError,match="non-negative"):
        Polyform.from_networkx(graph)


def test_parallel_sum_overflow_is_rejected():
    """Finite individual edges must not produce an infinite combined coefficient."""
    graph=nx.MultiGraph()
    graph.add_edge(0,1,weight=1e308)
    graph.add_edge(0,1,weight=1e308)
    with pytest.raises(ValueError,match="parallel"):
        Polyform.from_networkx(graph)


def test_mixed_and_colliding_node_labels_remain_distinct():
    """Integer 1, string '1' and a reserved suffix name must never collapse."""
    labels=[1,"1","1#1","",("tuple",2)]
    graph=nx.path_graph(labels)
    p,points=Polyform.from_networkx(graph,return_points="dict")
    assert list(points)==labels  # Mixed types use graph insertion order.
    assert len(set(points.values()))==len(labels)
    assert len({str(x) for x in points.values()})==len(labels)
    assert points[1] != points["1"]
    # Unit path P5 has endpoint degree one and interior degree two.
    expected=np.diag([1,2,2,2,1])-np.eye(5,k=1)-np.eye(5,k=-1)
    np.testing.assert_array_equal(p.to_matrix()[0],expected)


def test_comparable_labels_retain_sorted_order():
    """Preserve the established numeric ordering for existing graph callers."""
    graph=nx.Graph()
    graph.add_nodes_from([2,0,1])
    p,points=Polyform.from_networkx(graph,return_points="list")
    assert [str(x) for x in points]==["0","1","2"]
    assert p.vertices==tuple(points)


def test_network_api_rejects_invalid_graph_and_return_mode():
    """Bad constructor arguments produce explicit errors instead of silent fallbacks."""
    with pytest.raises(TypeError,match="NetworkX"):
        Polyform.from_networkx([])
    with pytest.raises(ValueError,match="return_points"):
        Polyform.from_networkx(nx.Graph(),return_points="mapping")


def test_signed_coefficients_remain_valid_outside_network_import():
    """Network validation must not forbid subtraction or signed algebraic expressions."""
    a,b=Point.create_list("ab")
    edge=quForm([(a,b)])
    signed=Polyform.from_pairs([([(a,b)],-2)])
    assert signed == -2*edge
    assert signed+3*edge == edge.as_polyform()
    assert signed.exp() == Polyform.one()-2*edge


def test_laplacian_clips_tiny_sign_noise():
    """A negative conductance within the explicit roundoff bound becomes zero."""
    delta=1e-12
    p=Polyform.from_laplacian([[-delta,delta],[delta,-delta]],Point.create_list("ab"))
    assert p.is_zero()
    assert len(p.vertices)==2


def test_laplacian_row_sum_check_has_no_relative_tolerance():
    """A wrong diagonal remains an error even when the conductances are large."""
    with pytest.raises(ValueError,match="row sums"):
        Polyform.from_laplacian([[1e8+1,-1e8],[-1e8,1e8+1]],Point.create_list("ab"))


def test_unrepresentably_large_weight_is_rejected():
    """A real integer too large for floating coefficients gets a clear validation error."""
    graph=nx.Graph()
    graph.add_edge(0,1,weight=10**1000)
    with pytest.raises(ValueError,match="finite"):
        Polyform.from_networkx(graph)
