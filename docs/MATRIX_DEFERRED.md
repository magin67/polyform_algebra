# Matrix module: deferred behavior and implementation stages

Updated 2026-10-30. Implemented behavior is specified in [MATRIX_API.md](MATRIX_API.md).
This ledger distinguishes mathematical decisions from scheduled implementation.

## Behavior requiring clarification before implementation

| Topic | Current boundary | Decision/evidence needed |
|---|---|---|
| Directed resistance and deviation | Basic symmetric distance formula is allowed only for vector-type matrices. On nonsymmetric input it loses the skew part. Directed semantic conversions raise NotImplementedError. | Choose symmetric resistance versus directional pair data, normalization and information to preserve; audit balanced/one-way notes. Do not infer invertibility for arbitrary directed data. |
| General directed Greenian | Directed Laplacians can be stored, classified and pseudoinverted numerically; the inverse remains generic. | Specify the inverse/gauge, left/right nullspaces and handling of reducible graphs. Balance is not silently imposed. |
| Degenerate metric projections | No projection API yet and no implicit pseudoinverse fallback. | Explicit quotient/gauge or selection rule where a projection is nonunique or nonexistent. |
| Metric-free orthonormal flag | Existing Basis.is_orthonormal is retained for old algebra APIs. Numerical kernels document their Euclidean coordinate metric. | Replace the old declaration only with a specified metric-aware API and migration of diagonal-form operations. |
| Space containment and affine structures | Space currently has identity and dimension; Basis owns it. | Represent embeddings, the affine-weight functional and affine origins without confusing a subspace with a subgraph. |
| Label and element mutation | Matrices detect basis-list membership/order changes. | Fully immutable algebra elements or versioned identity. Existing mutable Element names/weights must not be changed while used as dictionary keys. |

## Known mathematics, planned small implementation stages

| Capability | What remains | Acceptance example |
|---|---|---|
| Incidence objects | Explicit edge basis, oriented endpoints, conductances and adapter methods. Adjacency and its Laplacian conversions are implemented. | L=B diag(c) B.T and invariance under reversal of edge orientation. |
| Numerical coordinates | Dedicated bi-/di-coordinate results, Series adapters and explicit transition objects reusing Basis links. | Nonorthogonal basis: d=G b, with invertibility checked on the chosen span. |
| Metric restriction and unique projections | Embeddings, linear/affine subspaces and nonsingular restricted metric solves. | Orthogonal residual in the path midpoint example; singular restriction rejected. |
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
