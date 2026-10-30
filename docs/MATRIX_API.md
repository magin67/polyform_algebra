# Matrix API: implemented contracts

Updated 2026-10-30. This document supersedes the initial proposals in
[MATRIX_MODULE_DESIGN.md](MATRIX_MODULE_DESIGN.md). Unimplemented behavior and
open mathematical decisions are tracked in [MATRIX_DEFERRED.md](MATRIX_DEFERRED.md).

## One matrix class, explicit context

```python
from polyform_algebra import Basis, Space
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind

vertices = Space("vertices", dimension=3)
basis = Basis(["a", "b", "c"], space=vertices)
axis = MatrixAxis(basis=basis)
laplacian = LabeledMatrix(
    [[1, -1, 0], [-1, 3, -2], [0, -2, 2]],
    rows=axis, columns=axis,
    kind=MatrixKind.FORM_COEFFICIENTS,
    representation="laplacian",
)
greenian = laplacian.to_greenian()
resistance = laplacian.to_resistance()
recovered = resistance.to_laplacian()
```

Space owns a name and fixed dimension, but no metric. Space identity is object
identity: matching names do not identify independent spaces. Basis owns its
Space. MatrixAxis refers to that Basis and captures its order; its `space`
property delegates to Basis. It stores no second independent space value.

`Basis(elements)` remains available with unspecified space. Unique hashable
elements are now required in every constructor. This does not prove linear
independence of arbitrary algebra elements. `basis.bind_space(space)` attaches
an unspecified basis once. A basis of a fixed-dimensional Space cannot be
resized; create a new Space/Basis to describe another span. Legacy unbound bases
retain add/remove operations. `link_to()` accepts two bases of the same Space,
or two legacy unspecified bases; a specified and unspecified basis cannot be
silently linked. Existing frame-based column conventions are preserved.

The former `MatrixAxis(space="V", basis=basis)` call remains supported: it binds
the Basis to a Space. Reuse `basis.space` for other bases of that space. Two
independent calls with the same string do not establish an identification.
Direct mutation of the legacy `basis.elements` list invalidates matrix contexts.

Values are dense float64, copied on construction and backed by immutable bytes.
`to_numpy()` returns an independent copy. Operations return new matrices.
Generic matrices require finite data. Only validated resistance representations
allow positive infinity, and only between complete connected-component blocks.

## Role, representation and properties

`kind` retains the existing MatrixKind enum: OPERATOR, METRIC,
FORM_COEFFICIENTS, PAIRWISE. It records a role, not a proof of symmetry or
positive definiteness. A metric role can carry an indefinite matrix.

`representation` accepts the strings (or MatrixRepresentation members)
`generic`, `adjacency`, `bipartite`, `laplacian`, `greenian`, `resistance`. Meaning is never inferred from
values. `with_representation(name, kind=...)` explicitly declares and validates
meaning on a new object. These are representations of one class, not subclasses.

| Declaration | Validation and interpretation |
|---|---|
| adjacency | Finite nonnegative square conductances with exactly zero diagonal; directed data is allowed. Role is PAIRWISE. |
| bipartite | Finite nonnegative conductance table between two distinct Space instances, possibly rectangular. Role is PAIRWISE. |
| laplacian | Finite square outgoing Laplacian: zero row sums, nonpositive off-diagonal entries, nonnegative diagonal. Directed matrices are accepted. Role is OPERATOR or FORM_COEFFICIENTS. |
| greenian | Symmetric centered PSD matrix whose pseudoinverse is a nonnegative-conductance network Laplacian, with rank compatible with components. Role is METRIC or FORM_COEFFICIENTS. |
| resistance | Symmetric nonnegative zero-diagonal finite blocks, each reconstructing a connected network. Positive infinity separates blocks. Role is PAIRWISE. |
| generic | Shape and numerical storage only; no network claim. |

Named representations except bipartite require the same Space and ordered labels on both axes.
If labels are algebra Elements, they must have affine weight one. Empty networks
are allowed; singleton components have zero Laplacian/Greenian/self-resistance.

`vector_type(atol=None, rtol=None)` computes exactly these four cases in the
current coordinates, including rectangular matrices:

| 1.T A = 0 | A 1 = 0 | Result |
|---|---|---|
| yes | yes | vector |
| yes | no | left_vector |
| no | yes | right_vector |
| no | no | non_vector |

The zero matrix is vector type. This uses ones, not the affine weights of a
mixed basis. Affine-vector interpretation outside a point basis needs the
weight functional and is deliberately not inferred by this method.

Other numerical queries: `is_symmetric()`, `is_skew_symmetric()`, `rank()`,
`is_singular()` (square only), `kernel(side="right"|"left")`,
`is_positive_definite()`, `is_positive_semidefinite()`, `is_laplacian()`.
Definiteness checks require symmetry. Kernel columns are Euclidean-orthonormal
in the current coordinates, not orthonormal in an unspecified PMG metric.

## Basic and semantic transformations

`pseudoinverse()` computes the Moore-Penrose inverse using SVD in Euclidean
coordinates, exchanges axes and returns a generic OPERATOR. It is an explicit
numerical operation, not a projection policy or automatic network conversion.

`distance_transform()` requires matching axes and **vector type on both sides**.
It computes D(A)_ij=A_ii+A_jj-A_ij-A_ji and returns generic PAIRWISE data.
For symmetric A this is A_ii+A_jj-2A_ij. For a nonsymmetric vector matrix,
the skew part is lost; no directed resistance interpretation is asserted.

`deviation_transform()` requires finite symmetric zero-diagonal data and computes
-J R J/2, J=I-11.T/n, returning generic FORM_COEFFICIENTS. These transforms are
mutually inverse on symmetric vector matrices and symmetric hollow matrices.
Neither transform by itself implies positive definiteness or network realizability.

Semantic dispatch:

| Source | to_greenian() | to_resistance() | to_laplacian() |
|---|---|---|---|
| adjacency | Through Laplacian (symmetric only) | Through Laplacian (symmetric only) | diag(A1)-A, also directed |
| symmetric Laplacian | Componentwise pseudoinverse | Through Greenian | New equivalent object |
| Greenian | New equivalent object | Componentwise distance | Componentwise pseudoinverse |
| resistance | Componentwise deviation | New equivalent object | Through Greenian |
| generic | Explicit declaration required | Explicit declaration required | Explicit declaration required |

Every result is validated. In particular a centered PSD matrix can fail network
reconstruction because its inverse has a positive off-diagonal entry. A general
Euclidean distance table is not necessarily network resistance.

Semantic conversions return canonical roles: FORM_COEFFICIENTS for Laplacian
and Greenian, PAIRWISE for resistance and adjacency. Thus converting an operator-role
Laplacian to Laplacian preserves its numbers and vertex representation but
explicitly presents them as form coefficients.

Disconnected networks have one centroid gauge per component. `components`
returns tuples of vertex labels and retains isolates. Cross-component resistance
is +inf, never the finite result of applying D to unrelated Greenian blocks.
Semantic conversions process finite blocks. Basic numerical operations reject
infinite data; they cannot silently combine those blocks. For directed
Laplacians, conversion to themselves is supported; Greenian/resistance conversion
raises NotImplementedError until directed conventions are implemented.

`to_adjacency()` returns a new adjacency matrix from adjacency or any supported
Laplacian conversion. Off-diagonal entries are -L_ij and the diagonal is zero.
Self-loops are unsupported because L=diag(A1)-A cannot recover them. Negative
conductances arising from tolerated Laplacian sign errors are clipped to zero;
positive conductances, however small, are retained. Roundtrips are therefore
exact algebraically for valid loopless inputs, and numerical within tolerance
for inverse-derived Laplacians.

## Bipartite participation tables

`H.to_unipartite(space, representation="laplacian")` retains one axis Space
and eliminates the opposite part. Select the actual `H.rows.space` or
`H.columns.space` instance; equal names do not identify spaces. The input must
declare `representation="bipartite"`. Output may be `laplacian` (default) or
`adjacency`, with the retained axis on both sides in its original order.

For retained rows P and eliminated columns X, d_x=sum_p H_px and
A_pq=sum_x H_px H_qx/d_x for p != q. The diagonal of A is zero.
Equivalently L_P=diag(H1)-H diag(d)^-1 H.T on nonempty columns. Retaining
columns applies the same formula to H.T. The implementation builds L from
off-diagonal couplings to avoid subtractive cancellation in its diagonal.
Column scaling avoids overflow in intermediate products and degrees; results
outside finite float64 range still raise an error.

Zero columns contribute nothing; zero retained rows remain isolated. Empty
parts and a single retained vertex are supported. Disconnected inputs are
processed without coupling their components. This is elimination preserving
effective resistances, not induced-subgraph selection or common-neighbor
counts H H.T. Counts, hours, or other data become conductances only through
the application's modeling choice. Column normalization changes that model.

Original activity totals remain available from the unchanged H: sums of its
rows/columns. They are not the reduced graph's degrees and cannot generally
be recovered from the output. No activity metadata is silently attached to
the reduced matrix. Keep H when both statistics are needed.

See the executable [bipartite example](../examples/bipartite_couplings.py).
The existing strict pandas adapter accepts `representation="bipartite"` with
explicit row and column contexts; `to_dataframe()` exports the reduced matrix.
For NetworkX use `to_networkx(result.to_laplacian())`.

Source: [Couplings Within One Part of a Bipartite Graph](https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Couplings%20Within%20One%20Part%20of%20a%20Bipartite%20Graph.md),
equations (3), (6) and calculation limitations, v5 revision
`f4e80ba6d7bec8d3e32774b957fe921bf76de264`. Russian and English reduction
formulas and worked examples agree. Adjacency conventions follow the same
revision's Introduction note, “Laplacian: graph, electrical network and
quadratic form”, section “The Laplacian matrix”.

## Tolerances and conditioning

`MatrixTolerance(atol=1e-10, rtol=1e-12)` is immutable and stored on each matrix.
Queries and pseudoinverse may override the two values for that operation;
semantic conversions use the stored policy. Tolerances must be finite,
nonnegative real numbers, excluding booleans.

The cutoff is `atol + rtol * scale`. Row/column zero-sum tests use each row's
or column's absolute sum as reference scale; symmetry and sign checks use the
largest absolute entry; SVD rank, kernel and inverse share the largest singular
value as reference. Tolerances are dimensional and must be chosen for the data.
Small but genuine conductances are not deleted from the Laplacian's topology.
An unresolved component rank causes an explicit error rather than conversion
into disconnected vertices. Extremely ill-conditioned inverse pairs may require
a smaller tolerance or higher-precision work outside this float64 layer.

## Ordering and basis changes

`reordered(rows=..., columns=...)` requires full permutations. It preserves
Space identities, role and representation; named networks require the same
permutation on both axes, except bipartite tables whose two parts can be
permuted independently. Selecting a subset, taking a subgraph, and eliminating
vertices are separate operations, not aliases of reordering.

`change_basis(rows=axis, columns=axis, row_transition=M_r,
column_transition=M_c)` takes explicit invertible old-to-new coordinate maps,
consistent with Basis.link_to: x_new=M x_old. It applies:

- OPERATOR: M_r A M_c^-1.
- METRIC: M_r^-T G M_c^-1.
- FORM_COEFFICIENTS: M_r C M_c.T.
- PAIRWISE: rejected; only reordering is supported.

The target bases must belong to the original spaces. Arbitrary basis change
clears network representation, because point-coordinate sign/centering conditions
need not survive. Numerical coefficient operations do not extend formal polyform
algebra to unsupported bilinear forms.

## Interchange

`matrix.pandas.from_dataframe(..., rows, columns, kind, representation="generic",
tolerance=None)` requires exact matching labels/order; no automatic alignment.
`to_dataframe()` exports detached values, including +inf for disconnected
resistance. Context must be supplied on re-import; DataFrame metadata is not
trusted as a mathematical declaration. pandas stays optional and is imported
only inside adapter calls. Tuple IDs remain scalar Index labels.

`matrix.networkx.from_networkx(graph, vertices=None, weight="weight", space=None,
tolerance=None)` imports a simple graph or digraph as a Laplacian. A supplied
vertex list must be a full permutation; all isolated nodes survive. Existing
edges require finite nonnegative real weights, missing weights mean one and
zero weights mean no coupling. Self-loops and multigraphs are explicitly rejected.
Outgoing arc weights occupy adjacency rows; L=diag(A1)-A.

`matrix.networkx.to_networkx(laplacian)` exports Graph for symmetric matrices,
DiGraph otherwise. This describes the matrix, not the provenance of reciprocal
arcs. Every positive conductance and every vertex is retained. Conversion from
other representations is explicit through `to_laplacian()`.

Run `python -m examples.matrix_transformations` for the weighted path and an
isolated vertex. Numerical algorithms do not call polyform exponentials.

## Sources and validation scope

PMG `v5` revision **c27d4f781640c3639df8c0d6f1a9ca11795b8139**:

- [Mutual transformations](https://github.com/magin67/pmg-notes/blob/c27d4f781640c3639df8c0d6f1a9ca11795b8139/pmg-library/en/Introduction/Laplacian,%20Green%20Matrix,%20and%20Effective%20Resistance%20Matrix.md), transformation table, reconstruction conditions, weighted three-vertex path. Russian and English formulas agree.
- [Spaces and bases](https://github.com/magin67/pmg-notes/blob/c27d4f781640c3639df8c0d6f1a9ca11795b8139/pmg-library/en/Spaces/Spaces,%20Subspaces,%20and%20the%20Choice%20of%20a%20Basis.md), linear spaces, vector subspaces and bases.
- [Bilinear Laplacian](https://github.com/magin67/pmg-notes/blob/c27d4f781640c3639df8c0d6f1a9ca11795b8139/pmg-library/en/Balanced%20graphs/The%20Bilinear%20Laplacian%20of%20a%20Balanced%20Graph.md), outgoing convention, balance and directed 3-cycle. Numerical storage does not implement the formal bilinear extension in that note.

Tests verify inverse paths, Penrose equations, kernel equations, role-specific
invariance, componentwise resistances, preserved vertex order, explicit failures,
and legacy API regressions. They do not constitute proof of complete PMG coverage.
