"""Exact outgoing Laplacians and rooted forests for small directed graphs."""

from collections import Counter
from dataclasses import dataclass

from sympy import ImmutableMatrix, eye, ones, zeros

from .core.data import DomainError, InvalidInput, _labels, _number, _ordered


@dataclass(frozen=True, init=False)
class DirectedNetwork:
    """A loopless digraph with rational weights and ordered vertex labels."""

    labels: tuple
    arcs: tuple

    def __init__(self, labels, arcs):
        labels = _labels(labels)
        if not labels:
            raise InvalidInput("expected at least one vertex")
        index = {label: i for i, label in enumerate(labels)}
        weights = Counter()
        for arc in _ordered(arcs):
            arc = _ordered(arc)
            if len(arc) != 3:
                raise InvalidInput("arcs must be triples (source, target, weight)")
            source, target, weight = arc
            _labels((source, target))
            if source not in index or target not in index:
                raise InvalidInput("arc endpoints must be known vertices")
            weight = _number(weight, "exact")
            if weight < 0:
                raise InvalidInput("arc weights must be nonnegative")
            weights[index[source], index[target]] += weight
        object.__setattr__(self, "labels", labels)
        object.__setattr__(self, "arcs", tuple(
            (labels[i], labels[j], weight) for (i, j), weight in sorted(weights.items()) if weight
        ))

    def laplacian(self):
        """Return the exact matrix diag(outgoing weights) minus adjacency."""
        index = {label: i for i, label in enumerate(self.labels)}
        matrix = zeros(len(index))
        for a, b, weight in self.arcs:
            i, j = index[a], index[b]
            matrix[i, i] += weight
            matrix[i, j] -= weight
        return ImmutableMatrix(matrix)

    def tree_potentials(self):
        """Return the principal cofactors, in vertex order, for trees toward roots."""
        matrix = self.laplacian()
        return tuple(matrix.minor_submatrix(i, i).det() for i in range(matrix.rows))

    def balance(self):
        """Scale outgoing arcs by positive mean-one tree potentials."""
        potentials = self.tree_potentials()
        if any(value <= 0 for value in potentials):
            raise DomainError("positive balancing requires a strongly connected graph")
        mean = sum(potentials) / len(potentials)
        factors = tuple(value / mean for value in potentials)
        by_source = dict(zip(self.labels, factors))
        balanced = DirectedNetwork(
            self.labels, tuple((a, b, weight * by_source[a]) for a, b, weight in self.arcs)
        )
        return balanced, factors

    def centered_greenian(self):
        """Invert the Laplacian on zero-sum vectors of a connected balanced graph."""
        matrix = self.laplacian()
        n = matrix.rows
        uniform = ones(n) / n
        if ones(1, n) * matrix != zeros(1, n) or matrix.rank() != n - 1:
            raise DomainError("a connected balanced graph is required")
        return ImmutableMatrix((matrix + uniform).inv() - uniform)

    def one_way_resistances(self):
        """Return r_ij = Gr_jj - Gr_ij in the declared vertex order."""
        greenian = self.centered_greenian()
        return ImmutableMatrix(ones(greenian.rows, 1) * greenian.diagonal() - greenian)

    def vertex_gram(self):
        """Return the directional Gram matrix with local vertices and equal centroid pairings."""
        greenian = self.centered_greenian()
        diagonal = greenian.diagonal()
        offsets = ones(greenian.rows, 1) * diagonal
        return ImmutableMatrix(greenian.T - (offsets + offsets.T) / 2)

    def forest_toponome(self):
        """Enumerate forests directed toward roots, grouped by positive component type.

        Directions already determine each root, so the result is in vertex
        normalization. Enumeration is intended for small graphs.
        """
        positions = {label: i for i, label in enumerate(self.labels)}
        arcs = [(positions[a], positions[b], weight) for a, b, weight in self.arcs]
        counts = Counter()

        def collect(start, components, used_sources, weight):
            sizes = Counter(components).values()
            kind = tuple(sorted((size - 1 for size in sizes if size > 1), reverse=True))
            counts[kind] += weight
            for offset in range(start, len(arcs)):
                a, b, conductance = arcs[offset]
                first, second = components[a], components[b]
                if a not in used_sources and first != second:
                    merged = tuple(first if component == second else component for component in components)
                    collect(offset + 1, merged, used_sources | {a}, weight * conductance)

        collect(0, tuple(range(len(positions))), set(), 1)
        return dict(counts)
