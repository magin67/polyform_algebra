"""Topology shared by strict networks and the legacy numerical matrix API."""


def connected_components(labels, edges):
    """Return components in input order; every supplied edge is a real coupling.

    Callers validate weights first. This routine never uses a numerical cutoff
    and preserves isolated vertices, including heterogeneous immutable labels.
    """
    neighbors = {label: set() for label in labels}
    for u, v in edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    unseen, groups = set(labels), []
    for start in labels:
        if start not in unseen:
            continue
        found, pending = set(), [start]
        unseen.remove(start)
        while pending:
            u = pending.pop()
            found.add(u)
            for v in neighbors[u] & unseen:
                unseen.remove(v)
                pending.append(v)
        groups.append(tuple(x for x in labels if x in found))
    return tuple(groups)
