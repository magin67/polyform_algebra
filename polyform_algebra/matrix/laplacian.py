"""Validated Laplacians, resistance geometry and the undirected Dirichlet problem.

PMG Spaces, 'Selected Vertices and the Dirichlet Problem', equations (7)-(8):
solve L_II B = -L_IP. Rows of B are bi-coordinates of projections, not
expansions of the original internal vertices in the selected point basis.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from .._core.basis import Basis
from .._core.space import Space
from . import _numeric as num
from .context import LabeledMatrix, MatrixAxis, MatrixKind
from .results import (
    CircumsphereResult, SpanningTreeCentralityResult,
    EigenvectorCentralityResult, AnnihilatorsResult, BalancingResult, VertexValues,
)


@dataclass(frozen=True, eq=False, init=False)
class Laplacian(LabeledMatrix):
    """A LabeledMatrix with invariant laplacian representation.

    Directed outgoing Laplacians are supported as data. Each specialized
    operation checks its own domain; dirichlet currently requires symmetry.
    Non-Laplacian transformations return ordinary LabeledMatrix objects.
    """

    def __init__(self, values, *, rows, columns=None,
                 kind=MatrixKind.FORM_COEFFICIENTS,
                 representation="laplacian", tolerance=None):
        if representation != "laplacian":
            raise ValueError("Laplacian requires laplacian representation")
        super().__init__(values, rows=rows, columns=rows if columns is None else columns,
                         kind=kind, representation=representation, tolerance=tolerance)

    def annihilators(self, atol=None, rtol=None):
        """Return labeled numerical bases of ker(L.T) and ker(L).

        Under L=diag(outgoing degrees)-A, right annihilation means Lr=0,
        and left annihilation means l.T L=0. A single SVD and shared cutoff
        give both dimensions. Only a 1-D kernel has a sum-n representative;
        reducible digraphs may have zeros in their normalized left vector.
        This method neither balances the graph nor defines a directed metric.
        """
        a = self._array()
        n = len(a)
        tol = self._tol(atol, rtol)
        unit = num.scale(a)
        scaled = a / unit if unit else a.copy()
        u, singular_values, vh = np.linalg.svd(scaled, full_matrices=True)
        cutoff = tol.atol / unit + tol.rtol * num.scale(singular_values) if unit else 0.0
        rank = int(np.count_nonzero(singular_values > cutoff))
        if n and rank == n:
            raise ValueError("Laplacian kernel is unresolved at this tolerance")
        left_values, right_values = u[:, rank:].copy(), vh[rank:].T.copy()
        # Fix only a column's arbitrary sign. Rotations within multidimensional
        # nullspaces are noncanonical and may change with the numeric backend.
        for basis in (left_values, right_values):
            for j in range(basis.shape[1]):
                if basis[np.argmax(np.abs(basis[:, j])), j] < 0:
                    basis[:, j] *= -1
        floor = 64 * np.finfo(float).eps * max(1, n) * num.scale(singular_values)

        def check_residual(matrix, vectors):
            residual = matrix @ vectors
            norms = np.linalg.norm(residual, axis=0)
            limits = (cutoff + floor) * np.linalg.norm(vectors, axis=0)
            if not np.isfinite(residual).all() or np.any(norms > limits):
                raise ValueError("annihilator failed numerical residual checks")

        check_residual(scaled.T, left_values)
        check_residual(scaled, right_values)
        left = right = None
        if n - rank == 1:
            total = left_values[:, 0].sum()
            if abs(total) <= 64 * np.finfo(float).eps * max(1, n):
                raise ValueError("left annihilator sum normalization is numerically unresolved")
            normalized_left = left_values[:, 0] * (n / total)
            normalized_right = np.ones(n)
            check_residual(scaled.T, normalized_left[:, None])
            check_residual(scaled, normalized_right[:, None])
            left = VertexValues(normalized_left, axis=self.rows, name="left_annihilator")
            right = VertexValues(normalized_right, axis=self.columns, name="right_annihilator")

        def labeled_basis(values, side, vertex_axis):
            dimension = values.shape[1]
            axis = MatrixAxis(basis=Basis(
                [f"{side}_{i}" for i in range(dimension)],
                space=Space(f"{side} annihilator coefficients", dimension)))
            return LabeledMatrix(values, rows=vertex_axis, columns=axis,
                                 kind=MatrixKind.OPERATOR, tolerance=tol)

        return AnnihilatorsResult(
            labeled_basis(left_values, "left", self.rows),
            labeled_basis(right_values, "right", self.columns), left, right, tol)

    def balance(self, *, componentwise=False):
        """Balance outgoing weights by the positive, sum-n left annihilator.

        PMG 'Balancing a Directed Graph by Spanning-Tree Potentials', (1)-(3).
        Positive scaling preserves every arc and each source's allocation
        proportions. It changes the graph, not merely its coordinates.
        Each weak component must be strongly connected; componentwise mode
        returns separate results with independent mean-one normalizations.
        """
        a = self._array()
        if not isinstance(componentwise, bool):
            raise TypeError("componentwise must be a bool")
        off = a.copy()
        np.fill_diagonal(off, 0)
        # Generic Laplacian validation accepts roundoff-sized sign errors.
        # This operation promises support preservation, so do not clip them.
        if np.any(off > 0) or np.any(np.diag(a) < 0):
            raise ValueError("balancing requires nonnegative conductances without sign corrections")
        support = off < 0
        groups = num.components(support | support.T)
        if not componentwise and len(groups) != 1:
            raise ValueError("a nonempty connected graph is required; use componentwise=True")

        def reaches_all(arcs):
            reached, pending = {0}, [0]
            while pending:
                vertex = pending.pop()
                for neighbor in np.flatnonzero(arcs[vertex]):
                    if neighbor not in reached:
                        reached.add(neighbor)
                        pending.append(neighbor)
            return len(reached) == len(arcs)

        results = []
        for group in groups:
            arcs = support[np.ix_(group, group)]
            # Forward and reverse reachability from one vertex prove strong
            # connectivity. Never infer it from a nearly positive SVD vector.
            if not reaches_all(arcs) or not reaches_all(arcs.T):
                raise ValueError("positive support-preserving balancing requires strongly connected weak components")
            if len(groups) == 1:
                graph = self
            else:
                labels = tuple(self.rows.labels[i] for i in group)
                axis = MatrixAxis(basis=Basis(labels, space=Space(
                    f"{self.rows.space.name}: component {len(results)}", len(labels))))
                graph = Laplacian(a[np.ix_(group, group)], rows=axis,
                                  kind=self.kind, tolerance=self.tolerance)
            annihilators = graph.annihilators()
            if annihilators.left is None:
                raise ValueError("balancing annihilator is unresolved at this tolerance; choose a smaller tolerance")
            sigma = annihilators.left.to_numpy()
            if np.any(sigma <= 0):
                raise ValueError("positive balancing multipliers are numerically unresolved")
            # Identity on an already balanced input (within stored tolerance)
            # also prevents repeated balancing from introducing extra drift.
            if graph.vector_type() == "vector":
                sigma = np.ones(len(group))
            with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                values = sigma[:, None] * graph._array()
            if not np.isfinite(values).all() or not np.array_equal(values != 0, graph._array() != 0):
                raise ValueError("balancing exceeds numeric range or loses nonzero entries")
            balanced = Laplacian(values, rows=graph.rows, columns=graph.columns,
                                 kind=graph.kind, tolerance=graph.tolerance)
            if balanced.vector_type() != "vector":
                raise ValueError("balanced Laplacian failed zero row/column sum checks")
            results.append(BalancingResult(balanced, VertexValues(
                sigma, axis=graph.rows, name="balancing_multiplier")))
        return tuple(results) if componentwise else results[0]

    def _metric_components(self, componentwise):
        """Select finite connected metric domains without relating distinct components."""
        a = self._array()
        if not isinstance(componentwise, bool):
            raise TypeError("componentwise must be a bool")
        if not self.is_symmetric():
            raise NotImplementedError("this operation requires an undirected Laplacian")
        if not componentwise and len(self._components) != 1:
            raise ValueError("a nonempty connected graph is required; use componentwise=True")
        if len(self._components) == 1:
            return (self,)
        result = []
        for group in self._components:
            labels = tuple(self.rows.labels[i] for i in group)
            axis = MatrixAxis(basis=Basis(labels, space=Space(
                f"{self.rows.space.name}: component {len(result)}", len(labels))))
            result.append(Laplacian(a[np.ix_(group, group)], rows=axis,
                                    tolerance=self.tolerance, kind=self.kind))
        return tuple(result)

    def _tree_geometry(self):
        """Shared connected computation from a validated centered Greenian.

        PMG 'Centroid, Spanning-Tree, and Eigenvector Centralities', (4):
        delta_i=sum_j c_ij R_ij. The total degree is 2(n-1), not the sum of
        weighted graph degrees. Never multiply disconnected +inf by zero.
        """
        n = len(self.rows.labels)
        green = self.to_greenian().to_numpy()
        resistance = num.distance(green)
        adjacency = self.to_adjacency().to_numpy()
        degrees = (adjacency * resistance).sum(axis=1)
        eps = self.tolerance.threshold(max(1, n))
        if (not np.isfinite(degrees).all()
                or abs(degrees.sum() - 2 * (n - 1)) > eps
                or (n > 1 and (np.any(degrees < 1 - eps)
                               or np.any(degrees > np.count_nonzero(adjacency, axis=1) + eps)))):
            raise ValueError("expected tree degrees failed numerical consistency checks")
        return degrees, green, resistance

    def circumsphere(self, *, componentwise=False):
        """Return the resistance radius and labeled bi-coordinates of the center.

        Connected undirected graphs only by default. With componentwise=True,
        return a tuple, in component order, including singleton radius zero
        and center coordinate one. An empty graph then gives an empty tuple.
        This is the sphere in the affine hull, not a smallest enclosing ball
        or a new element of the isotropic extension.
        """
        graphs = self._metric_components(componentwise)
        results = []
        for graph in graphs:
            degrees, green, resistance = graph._tree_geometry()
            coordinates = 1 - degrees / 2
            r = np.diag(green)
            # rho = tr(Gr)/n + ||O-centroid||^2. This avoids constructing an
            # inverse of R and retains negative coordinates of the center.
            radius = float(r.mean() + (r @ graph._array() @ r) / 4)
            eps = graph.tolerance.threshold(max(1, len(r)))
            residual = resistance @ coordinates - 2 * radius
            metric_eps = graph.tolerance.threshold(max(num.scale(resistance), abs(2 * radius)))
            if (not np.isfinite(radius) or radius < 0
                    or abs(coordinates.sum() - 1) > eps
                    or not np.isfinite(residual).all()
                    or np.any(np.abs(residual) > metric_eps)):
                raise ValueError("circumsphere failed center or equal-radius numerical checks")
            results.append(CircumsphereResult(radius, VertexValues(
                coordinates, axis=graph.rows, name="center_coordinate")))
        return tuple(results) if componentwise else results[0]

    def spanning_tree_centrality(self, *, componentwise=False):
        """Return expected tree degrees and delta/[2(n-1)] normalized scores.

        Trees are sampled with probability proportional to the product of
        conductances. This is not a distribution over arbitrary forests.
        Componentwise results have separate normalizations; singleton scores
        are None because their trees contain no edge ends to sample.
        """
        graphs = self._metric_components(componentwise)
        results = []
        for graph in graphs:
            degrees, _, _ = graph._tree_geometry()
            n = len(degrees)
            scores = None if n == 1 else VertexValues(
                degrees / (2 * (n - 1)), axis=graph.rows, name="spanning_tree_centrality")
            results.append(SpanningTreeCentralityResult(VertexValues(
                degrees, axis=graph.rows, name="expected_tree_degree"), scores))
        return tuple(results) if componentwise else results[0]

    def eigenvector_centrality(self, *, componentwise=False):
        """Return the adjacency Perron eigenpair, with scores normalized by sum.

        For undirected connected graphs, use the largest algebraic eigenvalue,
        not the largest absolute eigenvalue (bipartite graphs also have -rho).
        Componentwise mode normalizes each component separately. A singleton
        has eigenvalue zero and score one; empty component mode returns ().
        """
        results = []
        for graph in self._metric_components(componentwise):
            n = len(graph.rows.labels)
            if n == 1:
                eigenvalue, scores = 0.0, np.ones(1)
            else:
                adjacency = graph.to_adjacency().to_numpy()
                unit = num.scale(adjacency)
                if unit == 0:
                    raise ValueError("connected adjacency has no numerically resolved coupling")
                # Scaling bounds the eigensolver input; restore the eigenvalue's
                # conductance units afterwards. Average only the asymmetry
                # already accepted by the matrix's symmetry tolerance.
                scaled = adjacency / unit
                symmetric = (scaled + scaled.T) / 2
                eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
                cutoff = graph.tolerance.atol / unit + graph.tolerance.rtol * num.scale(eigenvalues)
                if eigenvalues[-1] - eigenvalues[-2] <= cutoff:
                    raise ValueError("Perron eigenvalue gap is unresolved at this tolerance; choose a smaller tolerance")
                vector = eigenvectors[:, -1]
                # An eigenvector's global sign is arbitrary. Taking componentwise
                # absolute values could hide an incorrect or unresolved eigenpair.
                if vector.sum() < 0:
                    vector = -vector
                if not np.isfinite(vector).all() or np.any(vector <= 0):
                    raise ValueError("positive Perron vector is numerically unresolved")
                scores = vector / vector.sum()
                with np.errstate(over="ignore"):
                    eigenvalue = float(eigenvalues[-1] * unit)
                residual = scaled @ scores - eigenvalues[-1] * scores
                if (not np.isfinite(eigenvalue) or not np.isfinite(scores).all()
                        or np.any(np.abs(residual) > cutoff)):
                    raise ValueError("Perron eigenpair failed finite-value or residual checks")
            results.append(EigenvectorCentralityResult(eigenvalue, VertexValues(
                scores, axis=graph.rows, name="eigenvector_centrality")))
        return tuple(results) if componentwise else results[0]

    def eliminate(self, keep):
        """Eliminate internal vertices, preserving resistance on ordered keep.

        PMG 'Vertex Elimination and Induced Couplings', equations (4)-(6).
        Reuse the unique Dirichlet solve B=-solve(L_II, L_IP), then form
        L_PP+L_PI B. Every component must contain a retained vertex.
        The source is unchanged; directed elimination is not supported yet.
        """
        problem = self.dirichlet(keep)
        a = self._array()
        pi, ii = problem._selected_indices, problem._internal_indices
        if not ii:
            reduced = Laplacian(a[np.ix_(pi, pi)], rows=problem.influence_matrix.columns,
                                kind=self.kind, tolerance=self.tolerance)
            return EliminationResult(reduced, problem)
        with np.errstate(over="ignore", invalid="ignore"):
            response = a[np.ix_(pi, pi)] + a[np.ix_(pi, ii)] @ problem.influence_matrix._array()
        if not np.isfinite(response).all():
            raise ValueError("elimination response exceeds the finite numeric range")
        # Off-diagonal entries add induced conductances without cancellation.
        # Rebuild degrees from those couplings: direct diagonal subtraction
        # can leave a spurious response when only one vertex remains.
        values = response.copy()
        np.fill_diagonal(values, 0)
        values = values + (values.T - values) / 2
        with np.errstate(over="ignore", invalid="ignore"):
            np.fill_diagonal(values, -values.sum(axis=1))
        # Allow only corrections within the input row's numerical scale,
        # never conceal an inaccurate solve by repairing its row sums.
        scales = np.max(np.abs(a[list(pi)]), axis=1)
        if (not np.isfinite(values).all() or np.any(
                np.abs(values - response) > self.tolerance.threshold(scales)[:, None])):
            raise ValueError("elimination response failed symmetry or zero-sum checks")
        reduced = Laplacian(values, rows=problem.influence_matrix.columns,
                            kind=self.kind, tolerance=self.tolerance)
        return EliminationResult(reduced, problem)

    def dirichlet(self, selected):
        """Return influence and harmonic extension for an ordered vertex subset.

        Every connected component must meet the nonempty selected set.
        The internal principal block includes couplings to selected vertices;
        it is not an induced-subgraph Laplacian. Solve it directly, without
        a pseudoinverse or implicit choice of undetermined constants.
        """
        a = self._array()
        if not self.is_symmetric():
            raise NotImplementedError("Dirichlet currently requires an undirected Laplacian")
        if isinstance(selected, (str, bytes, set, frozenset)):
            raise TypeError("selected must be an ordered iterable of vertex labels")
        selected = tuple(selected)
        if not selected or len(set(selected)) != len(selected):
            raise ValueError("selected must be nonempty and contain distinct vertices")
        labels = self.rows.labels
        index = {label: i for i, label in enumerate(labels)}
        if any(label not in index for label in selected):
            raise ValueError("selected contains an unknown vertex")
        pi = tuple(index[label] for label in selected)
        chosen = set(pi)
        ii = tuple(i for i in range(len(labels)) if i not in chosen)
        internal = tuple(labels[i] for i in ii)
        for group in self._components:
            if not chosen.intersection(group):
                missing = tuple(labels[i] for i in group)
                raise ValueError(f"Dirichlet solution is not unique: component {missing!r} has no selected vertex")

        # Solve components independently, preserving exact zero influence
        # between them and avoiding artificial conditioning across scales.
        values = np.zeros((len(ii), len(pi)))
        output_rows = {i: k for k, i in enumerate(ii)}
        output_columns = {i: k for k, i in enumerate(pi)}
        for group in self._components:
            gi = [i for i in group if i not in chosen]
            gp = [i for i in pi if i in group]
            if not gi:
                continue
            block = a[np.ix_(gi, gi)]
            rhs = -a[np.ix_(gi, gp)]
            singular_values = np.linalg.svd(block, compute_uv=False)
            cutoff = self.tolerance.threshold(num.scale(singular_values))
            if np.any(singular_values <= cutoff):
                raise ValueError("Dirichlet internal block is unresolved at this tolerance; choose a smaller tolerance")
            solution = np.linalg.solve(block, rhs)
            eps = self.tolerance.threshold(1.0)
            if (not np.isfinite(solution).all() or np.any(solution < -eps)
                    or np.any(np.abs(solution.sum(axis=1) - 1) > eps)):
                raise ValueError("Dirichlet influence failed nonnegativity or unit-row-sum checks")
            # Do not clip or renormalize coefficients: that would alter the
            # solved equations and obscure conditioning or invalid input.
            values[np.ix_([output_rows[i] for i in gi],
                          [output_columns[i] for i in gp])] = solution

        selected_axis = MatrixAxis(basis=Basis(
            selected, space=Space(f"{self.rows.space.name}: selected vertices", len(pi))))
        internal_axis = MatrixAxis(basis=Basis(
            internal, space=Space(f"{self.rows.space.name}: internal vertices", len(ii))))
        influence = LabeledMatrix(values, rows=internal_axis, columns=selected_axis,
                                  kind=MatrixKind.OPERATOR, tolerance=self.tolerance)
        return DirichletResult(influence, self.rows, pi, ii)


@dataclass(frozen=True, eq=False)
class EliminationResult:
    """Response Laplacian and the Dirichlet data linking it to the source.

    dirichlet.influence_matrix contains bi-coordinates of internal projections;
    dirichlet.extend restores a harmonic field in the original vertex order.
    """

    laplacian: Laplacian
    dirichlet: "DirichletResult"


@dataclass(frozen=True, eq=False)
class DirichletResult:
    """Influence matrix and original-axis embedding of the selected/internal sets.

    Axis spaces are coordinate spaces of vertex subsets, not a new ambient
    metric or a declaration that internal vertices belong to the selected span.
    """

    influence_matrix: LabeledMatrix
    source_axis: MatrixAxis
    _selected_indices: tuple = field(repr=False)
    _internal_indices: tuple = field(repr=False)

    @property
    def selected(self):
        return self.influence_matrix.columns.labels

    @property
    def internal(self):
        return self.influence_matrix.rows.labels

    def extend(self, values):
        """Extend one scalar field; return a dict in the original vertex order.

        Accept an exact label-to-value mapping, or a 1-D list/tuple/ndarray
        in selected order. Convert a pandas Series explicitly with to_dict().
        This avoids discarding its labels or silently aligning missing values.
        """
        labels = self.source_axis.labels
        b = self.influence_matrix._array()
        if isinstance(values, Mapping):
            if set(values) != set(self.selected):
                raise ValueError("values must contain exactly the selected vertex labels")
            values = [values[label] for label in self.selected]
        elif not isinstance(values, (list, tuple, np.ndarray)):
            raise TypeError("values must be a mapping or positional list, tuple or ndarray")
        raw = np.asarray(values)
        if raw.dtype.kind not in "iuf":
            raise TypeError("prescribed values must be real numbers")
        f = np.array(raw, dtype=float)
        if f.shape != (len(self.selected),) or not np.isfinite(f).all():
            raise ValueError("prescribed values must be a finite vector in selected order")
        result = np.empty(len(labels))
        result[list(self._selected_indices)] = f
        result[list(self._internal_indices)] = b @ f
        if not np.isfinite(result).all():
            raise ValueError("harmonic extension exceeds the finite numeric range")
        return dict(zip(labels, result.tolist()))
