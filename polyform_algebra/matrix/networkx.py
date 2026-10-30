"""NetworkX adapters; numerical matrices retain the entire explicit vertex set."""

from numbers import Real
import numpy as np
from .._core.basis import Basis
from .._core.space import Space
from .context import MatrixAxis, MatrixKind, LabeledMatrix


def from_networkx(graph, *, vertices=None, weight="weight", space=None, tolerance=None):
    """Import a simple graph as an outgoing Laplacian, preserving isolated nodes.

    Directed graphs are accepted numerically. Missing weights mean one;
    zero weights create no coupling. Parallel-edge aggregation is deferred.
    """
    if graph.is_multigraph():
        raise NotImplementedError("parallel-edge aggregation is deferred; supply a simple graph")
    labels = list(graph.nodes) if vertices is None else list(vertices)
    if len(labels) != len(set(labels)) or set(labels) != set(graph.nodes):
        raise ValueError("vertices must be a full permutation, not an induced subgraph")
    n = len(labels)
    space = Space("network vertices", n) if space is None else space
    axis = MatrixAxis(basis=Basis(labels, space=space))
    index = {node: i for i, node in enumerate(labels)}
    adjacency = np.zeros((n, n))
    for u, v, data in graph.edges(data=True):
        if u == v:
            raise ValueError("self-loops are not supported")
        c = data.get(weight, 1) if weight is not None else 1
        if isinstance(c, bool) or not isinstance(c, Real) or not np.isfinite(c) or c < 0:
            raise ValueError("conductances must be finite nonnegative real numbers")
        adjacency[index[u], index[v]] = c
        if not graph.is_directed():
            adjacency[index[v], index[u]] = c
    lap = np.diag(adjacency.sum(axis=1)) - adjacency
    return LabeledMatrix(lap, rows=axis, columns=axis, kind=MatrixKind.FORM_COEFFICIENTS,
                         representation="laplacian", tolerance=tolerance)


def to_networkx(matrix, *, weight="weight"):
    """Export a Laplacian, preserving labels and all positive conductances.

    Symmetric data yields Graph, otherwise DiGraph. Matrix data alone cannot
    recover whether symmetric couplings originally came from reciprocal arcs.
    """
    import networkx as nx

    if matrix.representation != "laplacian":
        raise ValueError("convert explicitly to laplacian before graph export")
    a = matrix.to_numpy()
    undirected = matrix.is_symmetric()
    graph = nx.Graph() if undirected else nx.DiGraph()
    labels = matrix.rows.labels
    graph.add_nodes_from(labels)
    for i, u in enumerate(labels):
        for j, v in enumerate(labels):
            if i == j or (undirected and j <= i):
                continue
            c = -(a[i, j] + a[j, i]) / 2 if undirected else -a[i, j]
            if c > 0:
                graph.add_edge(u, v, **{weight: float(c)})
    return graph
