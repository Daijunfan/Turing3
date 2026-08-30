# Certificate-first contact-repair model

## Objects

- **Partial Scheme:** an explicit sum of exact rank-one products that need not
  yet equal its target.
- **Residual Completion Rank:** `kappa(T,P)=rank(T-P)`.
- **Contact Organ:** an exact scheme for a declared set of recursive slots,
  together with embeddings, field, commutativity semantics, baseline cost,
  and an independently verifiable certificate.
- **Contact-Rank Hypergraph:** slots are vertices; a positive-gain certified
  organ on slot set `Q` is a hyperedge.
- **Baseline Envelope `B(Q)`:** the minimum exact cost among independent
  evaluation, Kronecker composition, axis concatenation, serendipitous
  substitution, Pan pair fusion, projection/peeling, known direct-sum or
  aggregation constructions, and same-field public explicit certificates.
- **Excess Contact Gain:** `B(Q)-rank(organ)`. Only a strictly positive value
  creates a contact hyperedge.

Certificates have trust states `VERIFIED_EXPLICIT`, `VERIFIED_DERIVED`,
`METADATA_ONLY`, `INVALID`, `REPAIRED`, and `QUARANTINED`. Only the first two
enter compilation or search.

## Classification

Organs are classified as `INDEPENDENT`, `KRONECKER`, `AXIS_CONCAT`,
`SERENDIPITOUS`, `PAN_PAIR`, `DISJOINT_SUM`, `GAUGE_EQUIVALENT`,
`GENUINE_MIXED_CONTACT`, or `UNCLASSIFIED`. Failure to exclude a known
construction forces `UNCLASSIFIED`; it is not novelty evidence.

## Packing theorem

If certified hyperedges form an exact cover of the slot vertices and their
embeddings have disjoint declared slot ownership (or an explicitly certified
sharing rule), the union of their rank-one products computes the union target.

Proof: exactness holds on each covered slot set by its certificate. Disjoint
ownership prevents double counting. The exact cover exhausts every slot once,
so summing the local identities yields the parent identity. The
`ContactPackingSolver` checks disjointness and complete coverage before
returning a minimum-cost cover.

## Recursive complexity

A contact organ changes a recursion grammar edge multiplicity. For a
multi-type grammar, the exponent is the root of `rho(A(gamma))=1`. A rank
reduction improves the exponent only if the recomputed spectral root is
strictly smaller; a local rank saving by itself is insufficient.

## Relation to prior work

Flip graphs change rank-one decompositions while preserving their tensor.
Pan aggregation supplies a known cross-slot identity. Axis concatenation and
direct sums join compatible maps. Ordinary tensor decomposition ignores the
slot provenance used by the compiler. Tensor-rank subadditivity,
submultiplicativity, flattening bounds, Pan identities, and known ranks are
prior mathematics.

The project-specific contribution is the certificate-first workflow tying
source quarantine, residual localization, baseline-envelope comparison,
typed contact hyperedges, and exact-cover compilation together. This round
does not establish a new rank, identity, parameterized family, or exponent.
