# Matrix module: design and first foundation

**Current implementation (2026-10-30):** see [MATRIX_API.md](MATRIX_API.md) and
[deferred behavior](MATRIX_DEFERRED.md). The text below preserves the initial
2026-10-29 design and its first follow-up as historical rationale. Its proposed
classes, MatrixAxis-owned space, finite-only storage and unimplemented-operation
statements have been superseded by the current API: one matrix class, Basis-owned
Space, validated representations, componentwise conversions and explicit +inf
cross-component resistance. Do not use the historical tables as the current API.

Status: proposed architecture; only the storage/context and DataFrame adapters
below are implemented. Date: 2026-10-29. Code baseline: `d4c82cf` (254 tests).
Development branch: `codex/matrix-module`; merging into main is a separate decision.

## Decision: NumPy storage, explicit context, optional pandas

Use a plain dense real `numpy.ndarray` internally, wrapped by a small library
object. Use pandas DataFrame/Series at data boundaries, not as the mathematical
type system. NetworkX remains the graph interchange layer. None of the numerical
algorithms will require a polyform exponential.

DataFrame is useful for importing tables, attaching business identifiers and
exporting results. It cannot by itself distinguish a potential space from its
dual, an operator from a metric, or a list of vertices from an affine basis.
Its automatic alignment is useful in data analysis but must not silently
reorder a mathematical operation. NumPy structured arrays name record fields,
not both coordinate axes; they are not a replacement for context.

References: [pandas dot and alignment](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.dot.html),
[NumPy structured arrays](https://numpy.org/doc/stable/user/basics.rec.html),
[NetworkX NumPy conversion](https://networkx.org/documentation/stable/reference/generated/networkx.convert_matrix.to_numpy_array.html).
NetworkX conversion accepts an explicit `nodelist`; a partial list selects an
induced subgraph, so our future graph adapter must require a full permutation.

The initial storage uses float64, finite real entries and dense arrays. Exact
rationals, complex coefficients, sparse storage and out-of-core data are future
backend decisions. Large real datasets need sparse storage eventually; this
first dense interface makes no scalability claim. Missing values are rejected,
not replaced by zero. Unit conventions remain the caller's responsibility.

## Source audit

Read the introductory series and the eight-note Spaces series in `v5/pmg-library`.
Pinned source revision: `1e2a1221a4323b22e7fdbc96a4cb40254076cb5d`.
The following paths are relative to that revision, not the default branch:

- [Introduction index](https://github.com/magin67/pmg-notes/blob/1e2a1221a4323b22e7fdbc96a4cb40254076cb5d/pmg-library/en/Introduction/PMG%20Introductory%20Series.md): notes 1–4 define network matrices and resistance geometry; 5 and 9 Gram determinants; 6–7 variations; 8 and 10 separate formal algebra from metric evaluation.
- [Mutual transformations](https://github.com/magin67/pmg-notes/blob/1e2a1221a4323b22e7fdbc96a4cb40254076cb5d/pmg-library/en/Introduction/Laplacian,%20Green%20Matrix,%20and%20Effective%20Resistance%20Matrix.md): sections on centering, reconstruction and the weighted three-vertex path.
- [Spaces index](https://github.com/magin67/pmg-notes/blob/1e2a1221a4323b22e7fdbc96a4cb40254076cb5d/pmg-library/en/Spaces/Spaces%20and%20Subspaces%20in%20PMG.md): spaces, null vector, nonlocal elements, coordinates, contraction, Dirichlet problem, elimination and bipartite couplings.
- [Null vector](https://github.com/magin67/pmg-notes/blob/1e2a1221a4323b22e7fdbc96a4cb40254076cb5d/pmg-library/en/Spaces/The%20Null%20Vector,%20Isotropic%20Representation,%20and%20the%20Boundary.md): vertex and extended Gram matrices, locality, quotient by the null direction.
- [Coordinates](https://github.com/magin67/pmg-notes/blob/1e2a1221a4323b22e7fdbc96a4cb40254076cb5d/pmg-library/en/Spaces/Projections,%20Bi-coordinates,%20and%20Di-coordinates.md): di-coordinates, reciprocal bases and nondegenerate projections.
- [Elimination](https://github.com/magin67/pmg-notes/blob/1e2a1221a4323b22e7fdbc96a4cb40254076cb5d/pmg-library/en/Spaces/Vertex%20Elimination%20and%20Induced%20Couplings.md): response matrix, Schur complement, component conditions.
- [Bipartite couplings](https://github.com/magin67/pmg-notes/blob/1e2a1221a4323b22e7fdbc96a4cb40254076cb5d/pmg-library/en/Spaces/Couplings%20Within%20One%20Part%20of%20a%20Bipartite%20Graph.md): normalized contributions of shared neighbors and zero rows/columns.

Russian counterparts of the matrix transformation and Spaces notes were also
consulted. Translation metadata is not uniform: for example Spaces/Bases has
`translated_revision: 2` while Russian revision is 3; Coordinates records 8,
matching Russian revision 8. Do not infer synchronization from English revision
numbers. The inspected matrix/coordinate formulas agree; no result-changing
language conflict was identified. This is a focused interface audit, not a proof
of the full theory or a complete line-by-line translation comparison.

## Existing code and compatibility

`_core/basis.py: Basis` already supplies ordered elements, indices, frame-derived
transitions and inverse links. `Basis.link_to(other)` stores a matrix M whose
columns are the old elements in the other basis: x_other = M x_self. Preserve it.
Basis is mutable and does not validate independence or duplicates at construction.
Do not create a second basis class in this module.

`Polyform.from_networkx/from_laplacian` validate undirected nonnegative networks.
`Polyform.to_matrix()` returns a coefficient array and an element-to-index map.
Its support-derived vertex list cannot recover isolated vertices. It must keep
its return signature. `to_eigenbasis()` diagonalizes coefficients, not the
resistance embedding. `Polysimplex` already has coordinate/frame operations;
numerical adapters should use those, with explicit basis and grade checks.

Core must not import `matrix`. The matrix context can use Basis, and adapters
can use algebra objects. Importing the existing package still loads its existing
public imports; numerical independence means no exponential computation, not
that the old package initializer has become lazy.

### Basis centroid and affine-weight terminology

`Basis.centroid(name=None)` now returns `(element, frame)`, where the frame
is the arithmetic mean of the current basis elements with coefficients 1/n.
The returned element records that frame and its Basis; the basis is not changed.
This implementation supports nonempty bases of distinct grade-one `Element`
objects. Label-only numerical bases need an explicit adapter to algebra elements.

The weight is computed from the frame, epsilon(mean)=sum epsilon(h_i)/n.
For a point basis it is one, for a vector basis zero; for a mixed basis it
may be any real value. In the latter cases this method is an arithmetic-mean
extension of the centroid convention, not a claim that the result is a point.
No metric, sphere-center or isotropic-locality claim is made.

Use `Element(..., affine_weight=w)` and `element.affine_weight` for the value
of epsilon. Use `frame.affine_weight()` on a grade-one `Polysimplex`, including
the zero linear combination. This computes the scalar boundary, respecting
orientation. Point and Vector use weights one and zero respectively.

The old `Element(..., multiplicity=w)` keyword and `.multiplicity` attribute
share storage with the new name and remain available without warnings. Passing
both constructor names with different values raises ValueError. Existing
`LinearCombination.multiplicity()` also remains unchanged; it is legacy weighted
term metadata, not a general affine-weight functional. In particular, the
product metadata on higher-grade simplices and the -1/-2 zero/unit sentinels
must not be interpreted as affine weights. `affine_weight()` rejects scalars,
higher grades and formal forms. Fully removing this legacy storage requires a
separate internal refactor, not a global text substitution.

Terminology follows [PMG terminology, Affine weight and Earlier names](https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Foundations/PMG%20terminology.md)
and [Coordinates, Bi-coordinates and frames](https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Projections,%20Bi-coordinates,%20and%20Di-coordinates.md).
The weight functional on a mixed basis is weighted by epsilon(h_i); it is not
the sum of coordinates unless every basis element has weight one.

Subsequent design discussion also proposes moving space ownership into Basis
and deriving the space of MatrixAxis from its Basis. This has not been
implemented here. The metric-free `is_orthonormal` declaration likewise remains
a compatibility feature; a metric-aware replacement requires its own change.

## Implemented foundation

```python
from polyform_algebra import Basis
from polyform_algebra.matrix import MatrixAxis, MatrixKind, LabeledMatrix

vertices = MatrixAxis(space="example:vertex-values", basis=Basis(["a", "b", "c"]))
table = LabeledMatrix(
    [[1, -1, 0], [-1, 3, -2], [0, -2, 2]],
    rows=vertices, columns=vertices, kind=MatrixKind.FORM_COEFFICIENTS,
)
array = table.to_numpy()  # detached copy in the declared order
```

`MatrixAxis` captures a unique, hashable ordered sequence from the existing
Basis. It keeps that Basis and rejects subsequent changes to its element order
or membership when labels/data are requested. It does not certify independence.
Labels may be strings, tuples, numeric IDs or algebra Elements. Their hashes
and equality must remain stable. Empty axes are allowed. `space` is a nonempty
caller-controlled identifier; identical labels do not establish a mathematical
identification of spaces. No global registry or identity inference exists yet.

`MatrixKind` distinguishes OPERATOR, METRIC, FORM_COEFFICIENTS and PAIRWISE.
The generic container checks shape and finite real entries only. It does not
promise that a matrix marked METRIC is symmetric, or that a coefficient matrix
is a network Laplacian. No multiplication, projection or basis change is exposed
until role-specific checks exist. Reordering is also explicit and deferred.

`matrix.pandas.from_dataframe(frame, *, rows, columns, kind)` requires exact axis
labels and order; duplicates, missing values and nonnumeric values fail. A
caller who wants reordering must reindex the DataFrame explicitly first.
`to_dataframe(matrix)` returns independent storage. Tuple labels stay scalar
Index entries. Context cannot be reconstructed from exported values alone: the
import requires it explicitly. pandas is imported only when an adapter is used;
install it separately for now. A Series/coordinate adapter is deferred until
coordinates have a dedicated type. No dependency or packaging metadata changes.

## Proposed next interfaces (not implemented)

Use specialized checked objects over the foundation; representation names are
not interchangeable with MatrixKind. Before adding generic arithmetic, introduce
space/dual compatibility checks and shared numerical tolerance policy.

| Interface | Context and result contract |
|---|---|
| `NetworkMatrices.from_networkx(graph, *, vertices=None, weight="weight")` | Full explicit vertex permutation; isolated vertices preserved; reject directed graphs, loops, negative/nonfinite weights initially. Explicit policy for parallel edges; do not inherit an accidental conversion default. |
| `.adjacency()`, `.laplacian()`, `.incidence(*, oriented_edges)` | A symmetric, nonnegative, zero diagonal; L=diag(A1)−A. Incidence B has vertex rows, ordered oriented edge columns, −1 at tail, +1 at head. Store positive conductances c separately: L=B diag(c) B.T. Reversing an orientation changes B, not L. |
| `green(laplacian, *, components="raise")` | Connected default; Gr in centroid gauge, inverse on 1-perp and zero on constants. Explicit componentwise mode may return block Green data plus component IDs. |
| `resistance(green, *, components)` | R_ij=Gr_ii+Gr_jj−2Gr_ij within each component; cross-component resistances are infinite, never the finite values obtained by applying that formula to unrelated blocks. |
| `green_from_resistance(R)` / `laplacian_from_green(Gr)` | Validate symmetry, zero diagonal of R, centered PSD rank n−1, then inverse on 1-perp. Reconstructed L must have nonpositive off-diagonal entries. Euclidean distances alone do not guarantee a network. |
| `Coordinates(values, *, axis, kind="bi")` | Column coordinates; bi describes the actual object only if it belongs to the declared span. Di is the covector of inner products. Distinguish coordinates of a projection from coordinates of the original object. |
| `Transition(source, target, values)` | Explicit invertible coordinate map; adapter to Basis.transitions, same column convention. No mutation of frames just to import numerical data. |
| `MetricTensor(values, *, axis)` | Checked symmetric numerical metric, explicit signature/rank. Di=G bi; inverse conversion requires nonsingularity. A Gram matrix refers to an ordered family of vectors, which need not be independent. |
| `IsotropicModel.from_resistance(R, *, null_label)` | G=−R/2 on local vertices; Gz=[[G,1],[1.T,0]], ordered vertices then z. Distinguished z is nonzero, z²=0, a_i·z=1. Gr=J G J is a different Gram matrix. |
| `LinearSubspace(embedding)` / `AffineSubspace(origin, directions)` | Embedding E has independent direction columns in an ambient Basis. Restriction E.T G E. Affine origin and weight functional are explicit, especially in mixed bases. |
| `subspace.project(x, *, metric)` | Solve (E.T G E)b=E.T G x; fail explicitly if restricted metric is singular. Affine projection acts on x−origin. No silent pseudoinverse; a later explicit quotient/gauge interface may resolve nonuniqueness. |
| `network.induced_subgraph(keep)` | Select A_PP and recompute degrees, removing incident edges. |
| `network.eliminate(keep)` | L_PP−L_PI solve(L_II,L_IP), plus harmonic extension matrix and original ordering. Require a retained vertex in every participating component; uncovered components require an explicit separate policy. |
| `network.contract(partition)` | Full disjoint vertex partition, aggregation S, L'=S.T L S; combine conductances and discard internal loops. Distinct from elimination and metric restriction. |
| `bipartite_couplings(H, *, retained_axis, eliminated_axis)` | Rectangular nonnegative conductance table. For i≠j, c'_ij=sum_x H_ix H_jx/d_x, d_x=sum_i H_ix. Skip zero columns explicitly; keep zero rows and original degree totals. |

The bipartite result is the graph induced by *elimination* in the PMG sense,
not NetworkX's ordinary induced subgraph or an unweighted common-neighbor
projection. Its diagonal adjacency is zero, not the diagonal of H D^-1 H.T.
Preserve both input axes and original degrees for the two-part invariants.

For disconnected resistance output, the generic finite container cannot hold
infinity. Use a dedicated resistance result with component partition and masked
cross-component entries before exposing that mode; do not weaken finite storage
silently. Empty networks and singleton components need explicit conventions;
the singleton Laplacian/Green is zero, its self-resistance zero, tau=1.

## Basis-change conventions

All coordinate vectors are columns. Let M map old coordinates to new ones,
as Basis.link_to already does: x_new=M x_old.

- Operator V→U: A_new=M_U A_old M_V^-1.
- Metric evaluation x.T G y: G_new=M^-T G_old M^-1.
- Coefficients C in a sum of tensor products of basis elements:
  C_new=M C_old M.T.
- Bi-coordinates transform by M; di-coordinates transform by M^-T.

Thus a network's identical numerical L entries have two distinct uses: the
covariant energy matrix for potential arguments f.T L f, and coefficients of
the formal sum of vertex-difference squares. Their transformation laws differ
because the underlying arguments/spaces differ. The adapter must name the use.
Resistance pair tables have no arbitrary linear basis-change rule.

An isotropic nonzero vector does not imply a singular ambient metric: Gz is
nondegenerate for an affinely independent resistance simplex, with signature
(n,1). Its restriction to z-perp is degenerate with radical span(z); the quotient
is positive definite. Localizing a weight-one element by p−(p²/2)z is nonlinear
and distinct from orthogonal projection onto a linear or affine subspace.
Numerical bilinear evaluation does not add formal bilinear forms to the public
algebra. Polyform adapters may only return supported quadratic forms.

## Small implementation stages and acceptance checks

1. This change: numerical storage, Basis context, strict optional DataFrame
   interchange, example and validation tests. Existing API unchanged.
2. NetworkX adapter, A/B/L with orientation and ordering checks. Compare weighted
   chain and isolated vertices; no exponential dependencies. Then connected
   Gr/R, inverse reconstruction and explicit disconnected policy.
3. Coordinate/transition types and role-specific basis changes. Verify numerical
   invariance on nonorthogonal bases; use the existing Basis as the only basis
   implementation. Then metrics, Gram matrices and isotropic extension, including
   a nonsingular metric with a null vector and a singular restriction.
4. Linear/affine subspaces, restriction and unique projections; test orthogonal
   residuals and refusal of ambiguous projections. Verify path midpoint example.
5. Elimination/contraction and bipartite table reduction; compare Schur complement
   with star formulas, retained resistance and disconnected cases. Add pandas
   Series input only after coordinate semantics are fixed. Add sparse backend
   before claiming suitability for large graphs.

Every numerical mathematical test should cite its note/section and derive its
expected value. Storage tests exercise API safeguards rather than claiming PMG
theorems. Packaging was not changed in this first stage; when it changes, verify
wheel installation and tests outside the checkout as required by AGENTS.md.
