"""Compatibility regression for relocated adapters; no mathematical rules change."""
import importlib
import inspect

import networkx as nx
import numpy as np
import pytest


@pytest.mark.parametrize('old,new,names', [
    ('matrix.pandas','adapters.pandas',('from_dataframe','to_dataframe','to_series')),
    ('matrix.networkx','adapters.matrix_networkx',('from_networkx','to_networkx')),
])
def test_old_imports_are_identical_functions(old, new, names):
    """Compatibility paths expose the same implementation, signature and return contract."""
    a = importlib.import_module('polyform_algebra.'+old)
    b = importlib.import_module('polyform_algebra.'+new)
    for name in names:
        assert getattr(a,name) is getattr(b,name)
        assert inspect.signature(getattr(a,name)) == inspect.signature(getattr(b,name))
        assert getattr(a,name).__module__ == 'polyform_algebra.'+new


@pytest.mark.parametrize('directed', [False, True])
def test_matrix_graph_roundtrip_preserves_direction_order_and_isolate(directed):
    """Relocation must preserve outgoing arc weights and explicit isolated vertices."""
    from polyform_algebra.adapters.matrix_networkx import from_networkx, to_networkx
    graph = nx.DiGraph() if directed else nx.Graph()
    graph.add_nodes_from(['u','v','isolated'])
    graph.add_edge('u','v',weight=2)
    matrix = from_networkx(graph,vertices=['isolated','v','u'])
    restored = to_networkx(matrix)
    assert list(restored) == ['isolated','v','u']
    assert restored.is_directed() == directed
    assert restored['u']['v']['weight'] == 2
    assert restored.degree('isolated') == 0
    restored['u']['v']['weight'] = 99
    assert to_networkx(matrix)['u']['v']['weight'] == 2


def test_strict_and_matrix_graph_adapters_remain_distinct():
    """The 0.7 strict adapter keeps rational/parallel semantics instead of silently returning a matrix."""
    from fractions import Fraction as Q
    from polyform_algebra.core import Network, InvalidInput
    from polyform_algebra.adapters.networkx import from_networkx as strict
    from polyform_algebra.adapters.matrix_networkx import from_networkx as matrix
    graph = nx.MultiGraph()
    graph.add_edge('a','b',weight=Q(1,3))
    graph.add_edge('a','b',weight=Q(2,3))
    result = strict(graph,vertices=['a','b'])
    assert isinstance(result,Network) and result.edges == (('a','b',Q(1)),)
    with pytest.raises(NotImplementedError):
        matrix(graph)
    with pytest.raises(InvalidInput):
        strict(nx.DiGraph([(0,1)]),vertices=[0,1])


def test_relocated_dataframe_roundtrip_and_validation():
    """The table adapter retains axis identity, detached exports and strict label order."""
    pd = pytest.importorskip('pandas')
    from polyform_algebra import Basis, Space
    from polyform_algebra.matrix import MatrixAxis, MatrixKind
    from polyform_algebra.adapters.pandas import from_dataframe, to_dataframe
    rows = MatrixAxis(basis=Basis(['a','b'],space=Space('employees',2)))
    columns = MatrixAxis(basis=Basis(['x'],space=Space('documents',1)))
    frame = pd.DataFrame([[2],[1]],index=['a','b'],columns=['x'])
    matrix = from_dataframe(frame,rows=rows,columns=columns,
                            kind=MatrixKind.PAIRWISE,representation='bipartite')
    assert matrix.rows is rows and matrix.columns is columns
    exported = to_dataframe(matrix)
    np.testing.assert_array_equal(exported.to_numpy(),frame.to_numpy())
    frame.iloc[0,0] = 20
    exported.iloc[1,0] = 30
    np.testing.assert_array_equal(matrix.to_numpy(),[[2],[1]])
    with pytest.raises(ValueError):
        from_dataframe(frame.iloc[::-1],rows=rows,columns=columns,
                       kind=MatrixKind.PAIRWISE,representation='bipartite')
