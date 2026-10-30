"""PMG matrix identities and numerical/context regressions.

Source revision: c27d4f781640c3639df8c0d6f1a9ca11795b8139 (v5).
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Introduction/Laplacian,%20Green%20Matrix,%20and%20Effective%20Resistance%20Matrix.md
Sections 'Transformation table' and 'Example: a three-vertex path'.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/The%20Bilinear%20Laplacian%20of%20a%20Balanced%20Graph.md
Sections 'Conventions' and 'One directed 3-cycle'.
"""

import numpy as np
import pytest
from polyform_algebra import Basis, Space, Point
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind, MatrixTolerance


L = np.array([[1, -1, 0], [-1, 3, -2], [0, -2, 2.]])
GR = np.array([[1/2, -1/6, -1/3], [-1/6, 1/6, 0], [-1/3, 0, 1/3]])
R = np.array([[0, 1, 1.5], [1, 0, .5], [1.5, .5, 0]])


def matrix(a, representation="generic", kind=MatrixKind.FORM_COEFFICIENTS, tolerance=None):
    a = np.asarray(a)
    rows = MatrixAxis(basis=Basis(range(a.shape[0]), space=Space("V", a.shape[0])))
    cols = rows if a.shape[0] == a.shape[1] else MatrixAxis(
        basis=Basis(range(a.shape[1]), space=Space("U", a.shape[1])))
    return LabeledMatrix(a, rows=rows, columns=cols, kind=kind,
                         representation=representation, tolerance=tolerance)


@pytest.mark.parametrize("values, expected", [
    ([[1, -1], [-1, 1]], "vector"),
    ([[1, -1], [0, 0]], "right_vector"),
    ([[1, 0], [-1, 0]], "left_vector"),
    (np.eye(2), "non_vector"),
    (np.zeros((2, 3)), "vector"),
    (np.empty((0, 0)), "vector"),
])
def test_four_vector_types(values, expected):
    """Zero column/row sums correspond exactly to left/right annihilation by ones."""
    assert matrix(values).vector_type() == expected


def test_vector_classification_avoids_overflow():
    """Finite large positive entries cannot be vector type merely because sums overflow."""
    assert matrix(np.full((2, 2), 1e308)).vector_type() == "non_vector"


def test_directed_networkx_outgoing_convention():
    """A single arc a->b contributes [[c,-c],[0,0]], hence right_vector only."""
    import networkx as nx
    from polyform_algebra.matrix.networkx import from_networkx, to_networkx
    graph = nx.DiGraph()
    graph.add_edge("a", "b", weight=3)
    lap = from_networkx(graph)
    np.testing.assert_array_equal(lap.to_numpy(), [[3, -3], [0, 0]])
    assert lap.vector_type() == "right_vector"
    assert list(to_networkx(lap).edges) == [("a", "b")]
    with pytest.raises(ValueError, match="vector type"):
        lap.distance_transform()
    with pytest.raises(NotImplementedError):
        lap.to_resistance()


@pytest.mark.parametrize("bad", [-1, np.nan, np.inf, True, "2"])
def test_networkx_rejects_invalid_conductances(bad):
    """Nonnegative finite real conductances define the supported network input model."""
    import networkx as nx
    from polyform_algebra.matrix.networkx import from_networkx
    graph = nx.Graph()
    graph.add_edge(0, 1, weight=bad)
    with pytest.raises(ValueError, match="conductances"):
        from_networkx(graph)


def test_directed_cycle_distance_discards_circulation():
    """For the unit directed 3-cycle D has all off-diagonal entries 3; skew is lost."""
    a = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1.]])
    lap = matrix(a, "laplacian")
    assert lap.vector_type() == "vector"
    assert not lap.is_symmetric()
    distances = lap.distance_transform()
    np.testing.assert_allclose(distances.to_numpy(), 3 * (np.ones((3, 3)) - np.eye(3)))
    np.testing.assert_allclose(distances.deviation_transform().to_numpy(), (a + a.T)/2)
    assert distances.representation == "generic"
    with pytest.raises(NotImplementedError, match="directed"):
        lap.to_greenian()
    np.testing.assert_array_equal(lap.to_laplacian().to_numpy(), a)


@pytest.mark.parametrize("a", [[[1, -1], [0, 0]], [[1, 0], [-1, 0]], np.eye(2)])
def test_distance_requires_both_zero_sums(a):
    """User-defined domain excludes one-sided and non-vector matrices."""
    with pytest.raises(ValueError, match="vector type"):
        matrix(a).distance_transform()


def test_weighted_chain_all_semantic_paths():
    """Series conductances 1,2 give R12=1, R23=1/2, R13=3/2 (PMG example)."""
    lap = matrix(L, "laplacian")
    green = lap.to_greenian()
    resistance = lap.to_resistance()
    np.testing.assert_allclose(green.to_numpy(), GR, atol=1e-12)
    np.testing.assert_allclose(resistance.to_numpy(), R, atol=1e-12)
    np.testing.assert_allclose(green.to_laplacian().to_numpy(), L, atol=1e-12)
    np.testing.assert_allclose(resistance.to_greenian().to_numpy(), GR, atol=1e-12)
    np.testing.assert_allclose(resistance.to_laplacian().to_numpy(), L, atol=1e-12)
    np.testing.assert_allclose(green.to_resistance().to_numpy(), R, atol=1e-12)
    np.testing.assert_allclose(L @ green.to_numpy(), np.eye(3)-np.ones((3, 3))/3, atol=1e-12)
    assert green.rows is lap.rows
    assert green.representation == "greenian"
    assert resistance.kind == MatrixKind.PAIRWISE
    assert resistance.components == ((0, 1, 2),)


def test_abstract_distance_deviation_roundtrip():
    """Double centering inverts distance on symmetric vector matrices, even indefinite ones."""
    a = -L
    restored = matrix(a).distance_transform().deviation_transform()
    np.testing.assert_allclose(restored.to_numpy(), a, atol=1e-12)
    assert restored.representation == "generic"
    assert not restored.is_positive_semidefinite()


def test_disconnected_resistance_and_isolate_roundtrip():
    """A conductance-two edge has resistance 1/2; an isolated vertex has infinite cross-resistance."""
    a = np.array([[2, -2, 0], [-2, 2, 0], [0, 0, 0.]])
    lap = matrix(a, "laplacian")
    r = lap.to_resistance()
    assert r.components == ((0, 1), (2,))
    np.testing.assert_allclose(r.to_numpy(), [[0, .5, np.inf], [.5, 0, np.inf], [np.inf, np.inf, 0]])
    np.testing.assert_allclose(r.to_laplacian().to_numpy(), a, atol=1e-12)
    assert r.is_symmetric()
    for operation in (r.pseudoinverse, r.distance_transform, r.deviation_transform, r.vector_type):
        with pytest.raises(ValueError, match="finite"):
            operation()
    order = (2, 0, 1)
    permuted = r.reordered(rows=order, columns=order)
    assert permuted.components == ((2,), (0, 1))
    np.testing.assert_allclose(permuted.to_laplacian().to_numpy(), a[np.ix_(order, order)], atol=1e-12)


@pytest.mark.parametrize("n", [0, 1, 3])
def test_empty_and_edgeless_networks(n):
    """Singleton blocks have zero self-resistance; all distinct isolated vertices are disconnected."""
    lap = matrix(np.zeros((n, n)), "laplacian")
    r = lap.to_resistance().to_numpy()
    assert np.all(np.diag(r) == 0)
    assert np.all(np.isinf(r[~np.eye(n, dtype=bool)]))
    np.testing.assert_array_equal(lap.to_resistance().to_laplacian().to_numpy(), lap.to_numpy())


@pytest.mark.parametrize("values, rep, kind", [
    ([[1, 1], [1, 1]], "laplacian", MatrixKind.OPERATOR),
    ([[-1, 1], [1, -1]], "laplacian", MatrixKind.FORM_COEFFICIENTS),
    (np.eye(2), "greenian", MatrixKind.METRIC),
    ([[-1, 1], [1, -1]], "greenian", MatrixKind.METRIC),
    ([[0, 1, 5], [1, 0, 2], [5, 2, 0]], "resistance", MatrixKind.PAIRWISE),
    ([[0, 1, np.inf], [1, 0, 1], [np.inf, 1, 0]], "resistance", MatrixKind.PAIRWISE),
    ([[0, -np.inf], [-np.inf, 0]], "resistance", MatrixKind.PAIRWISE),
    (np.zeros((2, 2)), "resistance", MatrixKind.PAIRWISE),
    (R, "resistance", MatrixKind.METRIC),
])
def test_invalid_network_meaning_is_rejected(values, rep, kind):
    """Euclidean pair distances alone need not yield nonnegative network conductances."""
    with pytest.raises(ValueError):
        matrix(values, rep, kind)


def test_semantic_dispatch_requires_explicit_meaning():
    """The same symmetric vector array can have multiple meanings; no heuristic inference."""
    generic = matrix(L)
    assert generic.is_laplacian()
    with pytest.raises(ValueError, match="declare"):
        generic.to_greenian()
    np.testing.assert_allclose(generic.with_representation("laplacian").to_greenian().to_numpy(), GR, atol=1e-12)


def test_pseudoinverse_penrose_equations_and_kernels():
    """Rank-one rectangular A must satisfy all four Penrose equations and both kernel equations."""
    a = np.array([[1., 2, 3], [2, 4, 6]])
    m = matrix(a)
    inverse = m.pseudoinverse()
    p = inverse.to_numpy()
    np.testing.assert_allclose(a @ p @ a, a)
    np.testing.assert_allclose(p @ a @ p, p)
    np.testing.assert_allclose(a @ p, (a @ p).T)
    np.testing.assert_allclose(p @ a, (p @ a).T)
    np.testing.assert_allclose(a @ m.kernel(), 0, atol=1e-12)
    np.testing.assert_allclose(m.kernel("left").T @ a, 0, atol=1e-12)
    assert m.rank() == 1
    assert inverse.rows is m.columns
    assert inverse.representation == "generic"


def test_tolerance_cutoff_is_shared_and_does_not_silently_drop_edges():
    """A weak edge below the chosen cutoff must fail inversion, not become an isolate."""
    tiny = matrix([[1e-12, -1e-12], [-1e-12, 1e-12]], "laplacian")
    assert tiny.components == ((0, 1),)
    with pytest.raises(ValueError, match="unresolved"):
        tiny.to_greenian()
    a = matrix(np.diag([1., 1e-12]))
    assert a.rank() == 1
    assert a.rank(atol=0, rtol=0) == 2
    b = matrix([[1., -1 + 1e-9], [-1, 1]])
    assert b.vector_type(atol=1e-8) == "vector"
    assert b.vector_type(atol=1e-12, rtol=0) == "non_vector"


@pytest.mark.parametrize("value", [-1, np.nan, np.inf, True, "small"])
def test_invalid_tolerances(value):
    """All classification and factorization tolerances must be explicit finite nonnegative reals."""
    with pytest.raises(ValueError):
        MatrixTolerance(atol=value)
    with pytest.raises(ValueError):
        matrix(L).vector_type(rtol=value)


def test_basis_changes_preserve_operator_metric_and_contraction():
    """For x'=Mx, operator action, x.T G x and trace(G C) remain invariant."""
    g = matrix([[2., 1], [1, 3]], kind=MatrixKind.METRIC)
    c = LabeledMatrix([[4., 1], [1, 2]], rows=g.rows, columns=g.columns, kind=MatrixKind.FORM_COEFFICIENTS)
    a = LabeledMatrix([[1., 2], [3, 4]], rows=g.rows, columns=g.columns, kind=MatrixKind.OPERATOR)
    target = MatrixAxis(basis=Basis(["u", "v"], space=g.rows.space))
    m = np.array([[2., 1], [0, 1]])
    kwargs = dict(rows=target, columns=target, row_transition=m, column_transition=m)
    gn, cn, an = [q.change_basis(**kwargs).to_numpy() for q in (g, c, a)]
    x = np.array([3., -2])
    np.testing.assert_allclose((m @ x) @ gn @ (m @ x), x @ g.to_numpy() @ x)
    np.testing.assert_allclose(an @ (m @ x), m @ (a.to_numpy() @ x))
    np.testing.assert_allclose(np.trace(gn @ cn), np.trace(g.to_numpy() @ c.to_numpy()))


def test_space_identity_uniqueness_and_mutation_guards():
    """Space names are labels; distinct identities cannot be silently linked or substituted."""
    space = Space("V", 2)
    basis = Basis(["a", "b"], space=space)
    axis = MatrixAxis(basis=basis)
    assert axis.space is basis.space is space
    with pytest.raises(ValueError):
        Basis(["a", "a"])
    with pytest.raises(ValueError):
        Basis(["a"], space=space)
    with pytest.raises(ValueError):
        basis.bind_space(Space("V", 2))
    with pytest.raises(ValueError):
        basis.add_element("c")
    with pytest.raises(ValueError):
        basis.remove_element("a")
    with pytest.raises(ValueError):
        basis.link_to(Basis(["c", "d"], space=Space("V", 2)))
    assert basis.reordered(["b", "a"]).space is space


def test_context_mismatch_and_immutability():
    """Equal numbers and labels do not identify spaces; storage cannot be made writable."""
    m = matrix(L, "laplacian")
    other = MatrixAxis(basis=Basis(range(3), space=Space("V", 3)))
    with pytest.raises(ValueError, match="same space"):
        LabeledMatrix(L, rows=m.rows, columns=other, kind=m.kind, representation="laplacian")
    with pytest.raises(ValueError):
        m.reordered(rows=[2, 1, 0])
    with pytest.raises(ValueError):
        m._values.setflags(write=True)
    with pytest.raises(ValueError):
        m.rows.basis.reordered([0, 1])


def test_networkx_and_pandas_roundtrip_with_isolated_tuple_label():
    """Explicit order and isolated nodes must survive graph and DataFrame interchange."""
    import networkx as nx
    from polyform_algebra.matrix.networkx import from_networkx, to_networkx
    from polyform_algebra.matrix.pandas import from_dataframe, to_dataframe
    graph = nx.Graph()
    graph.add_nodes_from([("isolated", 0), "b", "a"])
    graph.add_edge("a", "b", weight=2)
    lap = from_networkx(graph, vertices=["b", ("isolated", 0), "a"])
    exported = to_networkx(lap)
    assert list(exported.nodes) == ["b", ("isolated", 0), "a"]
    assert exported["a"]["b"]["weight"] == 2
    with pytest.raises(ValueError):
        from_networkx(graph, vertices=["a", "b"])
    pytest.importorskip("pandas")
    r = lap.to_resistance()
    restored = from_dataframe(to_dataframe(r), rows=r.rows, columns=r.columns,
                              kind=r.kind, representation="resistance")
    np.testing.assert_allclose(restored.to_laplacian().to_numpy(), lap.to_numpy(), atol=1e-12)
