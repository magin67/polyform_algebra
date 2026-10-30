# Matrix module: deferred behavior and implementation stages

Updated 2026-10-30. Implemented behavior is specified in [MATRIX_API.md](MATRIX_API.md).
This ledger distinguishes mathematical decisions from scheduled implementation.

## Behavior requiring clarification before implementation

| Topic | Current boundary | Decision/evidence needed |
|---|---|---|
| Directed resistance and deviation | Basic symmetric distance formula is allowed only for vector-type matrices. On nonsymmetric input it loses the skew part. Directed semantic conversions raise NotImplementedError. | Choose symmetric resistance versus directional pair data, normalization and information to preserve; audit balanced/one-way notes. Do not infer invertibility for arbitrary directed data. |
| General directed Greenian | Directed Laplacians can be stored, classified and pseudoinverted numerically; the inverse remains generic. | Specify the inverse/gauge, left/right nullspaces and handling of reducible graphs. Balance is not silently imposed. |
| Degenerate metric projections | Undirected Dirichlet projections are implemented under uniqueness conditions; no general projection API or implicit pseudoinverse fallback. | Explicit quotient/gauge or selection rule where a general metric projection is nonunique or nonexistent. |
| Metric-free orthonormal flag | Existing Basis.is_orthonormal is retained for old algebra APIs. Numerical kernels document their Euclidean coordinate metric. | Replace the old declaration only with a specified metric-aware API and migration of diagonal-form operations. |
| Space containment and affine structures | Space currently has identity and dimension; Basis owns it. | Represent embeddings, the affine-weight functional and affine origins without confusing a subspace with a subgraph. |
| Label and element mutation | Matrices detect basis-list membership/order changes. | Fully immutable algebra elements or versioned identity. Existing mutable Element names/weights must not be changed while used as dictionary keys. |

## Known mathematics, planned small implementation stages

### Planned Laplacian methods

Laplacian and its undirected `dirichlet(selected)`, `circumsphere()` and
`spanning_tree_centrality()` methods are implemented. The latter two provide
explicit componentwise results; singleton centrality has scores=None.
The following methods are planned only; no placeholder calculations or
unimplemented public methods are introduced by these stages.

| Proposed method | Contract to implement | Acceptance checks |
|---|---|---|
| `annihilators()` | Labeled left/right nullspaces and their dimensions; expose a unique normalized vector only for a one-dimensional kernel. For outgoing L, L1=0; under strong connectivity normalize the positive left vector to sum n. Balance then gives ones on both sides. Reducible cases need explicit handling. | Lr=0 and l.T L=0; balanced examples; more than one component must not silently produce a preferred unique vector. |
| `eigenvector_centrality()` | Perron vector of adjacency reconstructed from L, normalized to sum one. Start with connected undirected graphs; define disconnected/singleton behavior explicitly. | A v=lambda_max v, sum v=1; the P5 vector proportional to (1,sqrt(3),2,sqrt(3),1). |

Spanning-tree centrality is not an alias for arbitrary forest centrality.
The latter requires a distribution over forests or a fixed component count
and a weight convention. Neither is silently inferred from a method name.

Sources for subsequent implementation:

- [Balancing a Directed Graph by Spanning-Tree Potentials](https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/Balancing%20a%20Directed%20Graph%20by%20Spanning-Tree%20Potentials.md).
- [Resistance Radius and Diameter of a Graph](https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Metric%20parameters/Resistance%20Radius%20and%20Diameter%20of%20a%20Graph.md).
- [Centroid, Spanning-Tree, and Eigenvector Centralities](https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Metric%20parameters/Centroid,%20Spanning-Tree,%20and%20Eigenvector%20Centralities.md).

### Further matrix capabilities

| Capability | What remains | Acceptance example |
|---|---|---|
| Incidence objects | Explicit edge basis, oriented endpoints, conductances and adapter methods. Adjacency and its Laplacian conversions are implemented. | L=B diag(c) B.T and invariance under reversal of edge orientation. |
| Numerical coordinates | General bi-/di-coordinate operations, Series import and explicit transition objects reusing Basis links. VertexValues and strict Series export already cover the implemented scalar vertex results. | Nonorthogonal basis: d=G b, with invertibility checked on the chosen span. |
| Metric restriction and unique projections | General embeddings, linear/affine subspaces and nonsingular restricted metric solves beyond the implemented graph Dirichlet problem. | Orthogonal residual in general coordinates; singular restriction rejected. |
| Isotropic extension | Distinguished z, vertex Gramian -R/2 and extended Gramian, locality and affine weight. | Nonsingular ambient metric with a null vector; degenerate restriction to z-perp. |
| Graph operations | Separate induced-subgraph, elimination and contraction methods. | Star elimination preserves leaf resistance; induced deletion and contraction produce distinct results. |
| Polyform adapters | Explicit Basis and vertex set for coefficient import/export; preserve isolates and signed Greenian coefficients. | Symmetric vector matrix roundtrip through supported pair-square forms. Directed bilinear export must remain unsupported. |
| Sparse storage | Select sparse backend and bounded-memory operations; dense inverses are not scalable. | Same small cases as dense backend plus memory checks on a large sparse graph. |
| Parallel graph edges | Explicit aggregation and provenance policy; currently rejected. | Sum of conductances with all node/edge ordering decisions documented. |
| Exact arithmetic | Optional rational/symbolic backend; current layer is float64. | Exact weighted-chain and near-singular reference results. |

These entries are not declarations that the corresponding mathematics is
unknown. They scope the next changes so that the entire original module roadmap
is not implemented in one patch. New implementations must read the current v5
note and move their contract into MATRIX_API.md with source-linked tests.
