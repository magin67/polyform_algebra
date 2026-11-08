"""Optional adapter checks; PMG v5 core definition: parallel conductances add."""
from fractions import Fraction as Q
import networkx as nx
import pytest
from polyform_algebra.core import InvalidInput
from polyform_algebra.adapters.networkx import from_networkx, to_networkx


def test_graph_adapter_exact_parallel_isolates_and_copy():
    """Two parallel weights 1/3 and 2/3 aggregate to 1 without losing an isolate."""
    graph = nx.MultiGraph()
    graph.add_nodes_from([1,'1','isolated'])
    graph.add_edge(1,'1',weight=Q(1,3))
    graph.add_edge(1,'1',weight=Q(2,3))
    net = from_networkx(graph,vertices=['isolated','1',1])
    assert net.edges == (('1',1,Q(1)),)
    copy = to_networkx(net)
    assert list(copy) == ['isolated','1',1]
    copy['1'][1]['weight'] = 5
    assert net.edges[0][2] == 1


@pytest.mark.parametrize('kind', ['directed','loop','missing','float'])
def test_graph_adapter_rejects_silent_conversion(kind):
    """Adapter domain violations must not coerce directed, looped or float exact input."""
    graph = nx.DiGraph() if kind == 'directed' else nx.Graph()
    graph.add_edge(0,0 if kind == 'loop' else 1,weight=0.5 if kind == 'float' else 1)
    with pytest.raises(InvalidInput):
        from_networkx(graph,vertices=[] if kind == 'missing' else list(graph))
