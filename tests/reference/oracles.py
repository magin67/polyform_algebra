"""Test-only references: edge-subset enumeration, grounded solves and determinants.

No polyform_algebra import is permitted here. These routines intentionally use
different representations from the production solvers and formal multiplication.
Sources and pinned examples are recorded in fixtures/reference_corpus_v1.json.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import factorial, prod

import sympy as sp


def forests(n, edges):
    """Enumerate labelled edge subsets; reject cycles by joining components."""
    for mask in range(1 << len(edges)):
        parent = list(range(n))
        chosen, weight = [], Q(1)

        def root(v):
            while parent[v] != v:
                v = parent[v]
            return v

        for i, (a, b, c) in enumerate(edges):
            if not mask & (1 << i):
                continue
            ra, rb = root(a), root(b)
            if ra == rb:
                break
            parent[ra] = rb
            chosen.append((a, b))
            weight *= Q(c)
        else:
            groups = defaultdict(list)
            for v in range(n):
                groups[root(v)].append(v)
            yield tuple(chosen), tuple(sorted(map(tuple, groups.values()))), weight


def forest_statistics(n, edges):
    """Rooted counts weight each forest by the product of component sizes."""
    ordinary, rooted = [Q(0)] * max(n, 1), [Q(0)] * max(n, 1)
    for chosen, groups, weight in forests(n, edges):
        k = len(chosen)
        ordinary[k] += weight
        rooted[k] += weight * prod(map(len, groups))
    return ordinary, rooted


def tree_and_separating_forests(n, edges):
    """R_ab = separating two-forest weight / tree weight, with no matrix solve."""
    tau = Q(0)
    separating = [[Q(0) for _ in range(n)] for _ in range(n)]
    for chosen, groups, weight in forests(n, edges):
        if len(groups) == 1:
            tau += weight
        elif len(groups) == 2:
            for a in groups[0]:
                for b in groups[1]:
                    separating[a][b] += weight
                    separating[b][a] += weight
    return tau, separating


def laplacian(n, edges):
    result = sp.zeros(n)
    for a, b, c in edges:
        c = sp.Rational(c)
        result[a, a] += c
        result[b, b] += c
        result[a, b] -= c
        result[b, a] -= c
    return result


def grounded_resistance(n, edges):
    """Ground the last vertex and solve all current columns; no pseudoinverse."""
    if n == 1:
        return sp.zeros(1)
    a = laplacian(n, edges)[:-1, :-1]
    h = a.LUsolve(sp.eye(n-1)).row_join(sp.zeros(n-1, 1)).col_join(sp.zeros(1, n))
    return sp.Matrix(n, n, lambda i, j: h[i, i]+h[j, j]-h[i, j]-h[j, i])


def cayley_menger(r, vertices):
    """Squared simplex volume from distances, independent of the production Gram path."""
    k = len(vertices)-1
    a = sp.ones(k+2)
    a[0, 0] = 0
    for i, u in enumerate(vertices):
        for j, v in enumerate(vertices):
            a[i+1, j+1] = r[u, v]
    return Q((-1)**(k+1)*a.det(), 2**k*factorial(k)**2)


def forest_tensor(n, edges, k):
    """Grade-k forest sum by incidence minors, without exterior multiplication."""
    indices = list(combinations(range(n), k))
    result = defaultdict(Q)
    for chosen in combinations(edges, k):
        columns = [sp.Matrix([int(i == b)-int(i == a) for i in range(n)]) for a, b, _ in chosen]
        matrix = sp.Matrix.hstack(*columns) if columns else sp.zeros(n, 0)
        coefficients = {I: Q(matrix.extract(I, range(k)).det()) for I in indices}
        weight = prod((Q(c) for _, _, c in chosen), start=Q(1))
        for I, x in coefficients.items():
            for J, y in coefficients.items():
                if x and y:
                    result[I, J] += weight*x*y
    return {key: value for key, value in result.items() if value}


def grove_weights(n, edges, selected):
    """Keep only forests whose every component meets the selected vertex set."""
    selected = set(selected)
    result = defaultdict(Q)
    for _, groups, weight in forests(n, edges):
        if all(selected.intersection(g) for g in groups):
            state = tuple(sorted(tuple(sorted(selected.intersection(g))) for g in groups))
            result[state] += weight
    return dict(result)


def minor_duality(a):
    """Jacobi complementary-minor identity, including orientation signs."""
    n, inverse = a.rows, a.inv()
    for k in range(n+1):
        for I in combinations(range(n), k):
            Ic = tuple(i for i in range(n) if i not in I)
            for J in combinations(range(n), k):
                Jc = tuple(j for j in range(n) if j not in J)
                assert a.extract(I, J).det() == (-1)**(sum(I)+sum(J))*a.det()*inverse.extract(Jc, Ic).det()
