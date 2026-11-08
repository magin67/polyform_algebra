"""Optional NetworkX exchange for strict networks, without floating conversion."""
from ..core import Network, InvalidInput, network


def from_networkx(graph, *, vertices, weight='weight', arithmetic='exact', policy=None):
    """Import an undirected graph with an explicit complete vertex order.

    Parallel edges are summed by the core; missing weights mean one.
    Loops, unsupported labels and negative conductances are rejected by PMG.
    """
    import networkx as nx
    if not isinstance(graph, nx.Graph) or graph.is_directed():
        raise InvalidInput('expected an undirected NetworkX graph')
    if not isinstance(vertices, (tuple, list)) or set(vertices) != set(graph.nodes):
        raise InvalidInput('vertices must explicitly include every node')
    return network(vertices, tuple((u, v, data.get(weight, 1))
                                  for u, v, data in graph.edges(data=True)),
                   arithmetic=arithmetic, policy=policy)


def to_networkx(net, *, weight='weight'):
    """Export a detached graph, preserving isolates, labels and rational weights."""
    import networkx as nx
    if not isinstance(net, Network):
        raise InvalidInput('expected Network')
    graph = nx.Graph()
    graph.add_nodes_from(net.context.labels)
    for u, v, c in net.edges:
        graph.add_edge(u, v, **{weight: c})
    return graph
