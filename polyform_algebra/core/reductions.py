"""Dirichlet extension and graph reductions for finite undirected networks.

Source: PMG v5 c909ddf8a590c69ed0a82fe48ca496e88fb707e0,
Resistance Geometry and Duality / Spaces: Selected Vertices and the Dirichlet
Problem (matrix notation); Vertex Elimination and Induced Couplings (response);
Vertex Contraction and Graph Gluing (vertex identification).
See docs/CORE_REDUCTIONS_API.md for links and conventions.
"""
from collections.abc import Mapping
from dataclasses import dataclass, field
from math import fsum, isfinite

import numpy as np
import sympy as sp

from .data import (
    Context, Coordinates, MatrixData, Network, InvalidInput, DomainError,
    NonuniqueProblem, NumericalUnresolved, _labels, _compatible,
)
from .network_math import SolveCheck, _require, _result, _finite, laplacian


@dataclass(frozen=True)
class DirichletResult:
    P: tuple
    I: tuple
    B: MatrixData
    H: MatrixData
    _source: Network = field(repr=False)
    _checks: tuple = field(default=(), repr=False)

    def extend(self, boundary_values):
        """Apply this extension; return Calculation with source-space coordinates."""
        return extend(self, boundary_values)


@dataclass(frozen=True)
class ReductionResult:
    network: Network
    H: MatrixData
    P: tuple
    I: tuple


@dataclass(frozen=True)
class ContractionResult:
    network: Network
    S: MatrixData


def _selection(net, labels):
    _require(net)
    labels = _labels(labels)
    if not set(labels).issubset(net.context.labels):
        raise InvalidInput('selected vertices must belong to the network')
    return labels


def _operator(values, rows, columns):
    return MatrixData(tuple(tuple(row) for row in values), rows, columns, 'operator')


def _solve(a, c, labels, policy):
    """Solve a grounded positive-definite block, with a scale-free error budget."""
    unit = float(np.max(np.abs(a)))
    if not unit or not isfinite(unit):
        raise NumericalUnresolved('grounded block has no finite positive scale')
    aa, cc = a / unit, c / unit
    if np.any((a != 0) & (aa == 0)) or np.any((c != 0) & (cc == 0)):
        raise NumericalUnresolved('grounded scaling would erase a nonzero coefficient')
    n = len(a)
    budget = max(policy.rtol, 100 * n * np.finfo(float).eps)
    try:
        ev = np.linalg.eigvalsh(aa)
        condition = float(ev[-1] / ev[0]) if ev[0] > 0 else float('inf')
        if ev[0] <= budget * ev[-1] or 10 * n * np.finfo(float).eps * condition > budget:
            raise NumericalUnresolved('grounded block is too ill-conditioned to certify')
        b = _finite(np.linalg.solve(aa, cc))
    except np.linalg.LinAlgError as error:
        raise NumericalUnresolved('grounded solve did not converge') from error
    denominator = np.linalg.norm(aa, np.inf) * np.linalg.norm(b, np.inf) + np.linalg.norm(cc, np.inf)
    residual = float(np.linalg.norm(aa @ b - cc, np.inf) / denominator)
    if residual > budget or np.any(b < 0) or np.max(np.abs(b.sum(axis=1) - 1)) > budget:
        raise NumericalUnresolved('harmonic residual or maximum principle could not be certified')
    # B_xp is positive exactly when p is reachable from x without passing
    # through another selected vertex. A zero caused by underflow would erase
    # an induced edge; conversely, roundoff must not invent a forbidden path.
    unseen = set(range(n))
    while unseen:
        pending = [next(iter(unseen))]
        group = set(pending)
        unseen.difference_update(group)
        while pending:
            x = pending.pop()
            reached = {int(y) for y in np.flatnonzero(a[x] < 0)} & unseen
            unseen.difference_update(reached)
            group.update(reached)
            pending.extend(reached)
        rows = list(group)
        reachable = np.any(c[rows] > 0, axis=0)
        if np.any((b[rows] > 0) != reachable):
            raise NumericalUnresolved('harmonic influence support was lost or invented numerically')
    return b, SolveCheck(labels, unit, n, condition, residual)


def dirichlet(net, boundary):
    """D01: B=-solve(L_II,L_IP), H_P=identity; each component needs an anchor.

    P follows the supplied order; I follows the source order. New P/I contexts
    are independent spaces, while H.rows retains the original context.
    No pseudoinverse chooses an arbitrary constant on unanchored components.
    """
    p = _selection(net, boundary)
    labels = net.context.labels
    if not labels:
        raise DomainError('Dirichlet extension is undefined for an empty network')
    if not p:
        raise InvalidInput('select at least one vertex')
    missing = tuple(group for group in net.components if not set(group).intersection(p))
    if missing:
        raise NonuniqueProblem('every component must contain a selected vertex', parameters=missing)
    internal = tuple(x for x in labels if x not in p)
    pc, ic = Context(p, net.context.arithmetic), Context(internal, net.context.arithmetic)
    idx = {x: i for i, x in enumerate(labels)}
    pi, ii = {x: i for i, x in enumerate(p)}, {x: i for i, x in enumerate(internal)}
    b = [[0] * len(p) for _ in internal]
    checks = []
    if internal:
        l = laplacian(net).value.to_numpy()
        # Solve components separately: unrelated large conductances must not
        # change the numerical rank of a small but well-conditioned component.
        for group in net.components:
            local_i = tuple(x for x in internal if x in group)
            local_p = tuple(x for x in p if x in group)
            if not local_i:
                continue
            ai, ap = [idx[x] for x in local_i], [idx[x] for x in local_p]
            a, c = l[np.ix_(ai, ai)], -l[np.ix_(ai, ap)]
            if net.context.arithmetic == 'exact':
                solved = sp.Matrix(a.tolist()).LUsolve(sp.Matrix(c.tolist()))
            else:
                solved, check = _solve(a, c, local_i, net.policy)
                checks.append(check)
            for i, x in enumerate(local_i):
                for j, y in enumerate(local_p):
                    b[ii[x]][pi[y]] = solved[i, j]
    h = [([int(x == y) for y in p] if x in pi else b[ii[x]]) for x in labels]
    value = DirichletResult(p, internal, _operator(b, ic, pc),
                            _operator(h, net.context, pc), net, tuple(checks))
    return _result(net, 'dirichlet', value, checks)


def extend(result, boundary_values):
    """D02: apply H to an exact label mapping or coordinates on H.columns."""
    if not isinstance(result, DirichletResult):
        raise InvalidInput('expected DirichletResult (the value of dirichlet())')
    if isinstance(boundary_values, Coordinates):
        _compatible(result.H.columns, boundary_values.context)
        values = boundary_values.values
    elif isinstance(boundary_values, Mapping):
        values = Coordinates(result.H.columns, boundary_values).values
    else:
        raise InvalidInput('boundary values require a label mapping or compatible Coordinates')
    numeric = result.H.rows.arithmetic == 'numeric'
    output = []
    for row in result.H.values:
        terms = [a * b for a, b in zip(row, values)]
        try:
            value = fsum(terms) if numeric else sum(terms)
        except (OverflowError, ValueError) as error:
            raise NumericalUnresolved('harmonic extension overflowed') from error
        if numeric and (not isfinite(value) or any(a != 0 and b != 0 and t == 0 for a, b, t in zip(row, values, terms))):
            raise NumericalUnresolved('harmonic extension lost a nonzero product or overflowed')
        output.append(value)
    return _result(result._source, 'extend', Coordinates(result.H.rows, tuple(output)), result._checks)


def _sum(weights, numeric):
    try:
        value = fsum(weights) if numeric else sum(weights)
    except OverflowError as error:
        raise NumericalUnresolved('conductance sum overflowed') from error
    if numeric and not isfinite(value):
        raise NumericalUnresolved('conductance sum overflowed')
    return value


def kron_reduce(net, keep):
    """D03: retain the response and pair resistances, using induced couplings.

    Off-diagonal conductances are c_pq + sum_x c_px B_xq. This avoids
    subtracting almost equal diagonal terms in the Schur complement.
    """
    calculation = dirichlet(net, keep)
    d = calculation.value
    neighbors = {x: {} for x in net.context.labels}
    for u, v, c in net.edges:
        neighbors[u][v] = neighbors[v][u] = c
    numeric = net.context.arithmetic == 'numeric'
    ii = {x: i for i, x in enumerate(d.I)}
    weights = []
    corrections = []
    for i, p in enumerate(d.P):
        for j in range(i + 1, len(d.P)):
            q = d.P[j]
            directions = []
            for u, v, k in ((p, q, j), (q, p, i)):
                terms = [neighbors[u].get(v, 0)]
                for x, c in neighbors[u].items():
                    if x in ii:
                        influence = d.B.values[ii[x]][k]
                        term = c * influence
                        if numeric and c and influence and not term:
                            raise NumericalUnresolved('induced conductance underflowed')
                        terms.append(term)
                directions.append(_sum(terms, numeric))
            a, b = directions
            if numeric:
                budget = max(net.policy.rtol, 100 * len(net.context.labels) * np.finfo(float).eps)
                if abs(a - b) > net.policy.conductance_atol + budget * max(a, b):
                    raise NumericalUnresolved('Kron response is not symmetric within policy')
                # The mean lies between the two responses, without a+b overflow.
                c = min(a, b) + abs(a - b) / 2
                if a != b:
                    corrections.append(('symmetric conductance', p, q, a, b, c))
            else:
                if a != b:
                    raise DomainError('exact Kron response is not symmetric')
                c = a
            if c:
                weights.append((p, q, c))
    reduced = Network(d.H.columns, tuple(weights), net.policy)
    # Rebuilding degrees from off-diagonal couplings enforces zero row sums.
    # Check the full response too: a small solve residual alone is insufficient.
    if numeric:
        l = laplacian(net).value.to_numpy()
        h = d.H.to_numpy()
        response = _finite(l @ h)
        expected = np.zeros_like(response)
        positions = [net.context.labels.index(x) for x in d.P]
        expected[positions, :] = laplacian(reduced).value.to_numpy()
        for group in net.components:
            rows = [net.context.labels.index(x) for x in group]
            scale = max(float(np.max(np.abs(l[rows]))), np.finfo(float).tiny)
            budget = max(net.policy.rtol, 100 * len(rows) * np.finfo(float).eps)
            if np.max(np.abs(response[rows] / scale - expected[rows] / scale)) > budget:
                raise NumericalUnresolved('Kron response residual exceeds policy')
    value = ReductionResult(reduced, d.H, d.P, d.I)
    return _result(net, 'kron_reduce', value, d._checks, corrections)


def contract(net, partition):
    """D04: identify a full ordered partition; S pulls new potentials back.

    Sum only cross-group edges. Huge internal conductances must never cancel
    small external ones through floating-point evaluation of S.T @ L @ S.
    """
    _require(net)
    if not isinstance(partition, Mapping):
        raise InvalidInput('partition must map new labels to ordered member sequences')
    groups = _labels(tuple(partition))
    owner = {}
    for group in groups:
        members = _labels(partition[group])
        if not members:
            raise InvalidInput('partition groups must be nonempty')
        for x in members:
            if x not in net.context.labels or x in owner:
                raise InvalidInput('partition contains an unknown or repeated vertex')
            owner[x] = group
    if set(owner) != set(net.context.labels):
        raise InvalidInput('partition must cover every vertex')
    context = Context(groups, net.context.arithmetic)
    index = {x: i for i, x in enumerate(groups)}
    buckets = {}
    for u, v, c in net.edges:
        a, b = owner[u], owner[v]
        if a != b:
            pair = tuple(sorted((a, b), key=index.__getitem__))
            buckets.setdefault(pair, []).append(c)
    edges = tuple((a, b, _sum(weights, net.context.arithmetic == 'numeric'))
                  for (a, b), weights in buckets.items())
    reduced = Network(context, edges, net.policy)
    s = _operator([[int(owner[x] == g) for g in groups] for x in net.context.labels], net.context, context)
    return _result(net, 'contract', ContractionResult(reduced, s))


def induced_subgraph(net, keep):
    """D05: retain only edges with both ends in keep; recompute vertex degrees."""
    keep = _selection(net, keep)
    context = Context(keep, net.context.arithmetic)
    edges = tuple((u, v, c) for u, v, c in net.edges if u in keep and v in keep)
    return _result(net, 'induced_subgraph', Network(context, edges, net.policy))
