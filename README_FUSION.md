# GermSynth-F

Proof-carrying cross-slot fusion extraction built on the GermSynth-R prototype.

## Reproduce

```bash
./reproduce_fusion.sh
```

The command reruns the GermSynth-R baseline, parses the pinned FMM-Lille Maple
Triads with exact rational coefficients, verifies every tensor coordinate,
recovers the outer/inner slot structure, constructs the product-slot
hypergraph, computes exact component residuals, generates JSON certificates,
and runs two standalone C++17 verifiers plus coefficient-mutation tests.

## Main result

The downloaded `8x27x30_tensor.mpl.bz2` and LRP declare 3736 products, while
the current FMM-Lille web page declares rank 3744 and the raw evaluation file
contains 3825 products. The explicit 3736 artifact is not an exact
matrix-multiplication scheme: it has 28,098 nonzero residual
coordinates modulo 1,000,003. A rational counterexample is
`A[1,25] * B[25,10] -> C[1,10]`, whose expected coefficient is 1 and whose
actual coefficient is 0.

Automatic structural analysis nevertheless recovers 238 ordinary slots and
six exceptional two-slot candidates. Five use 29 products and one uses 21,
so their nominal arithmetic is

```text
238*15 + 5*29 + 1*21 = 3736
250*15 - 3736 = 14
```

Every exceptional component has nonzero exact residual. Consequently no
rank-29 fusion kernel, valid gain-14 explanation, reusable constructor, new
rank, or new exponent is claimed.

See [results/HANDOFF_TO_BRAIN.md](results/HANDOFF_TO_BRAIN.md) for the complete
status table and [external/SOURCES.json](external/SOURCES.json) for pinned
commits, licenses, dates, and source-file hashes.

## Exact scheme format

Python uses sparse rows `{coordinate: Fraction}` and exports certificate rows
as `[coordinate, numerator, denominator]`. The unified axes are row-major
`A[m,n]`, `B[n,p]`, and `C[m,p]`. Importers normalize FMM-Lille's explicit
`p x m` output matrices and Perminov's `C^T` convention.

## No-contact proposition

Within ordinary non-overlap composition, suppose every rank-one product is
incident to at most one theoretical slot. The product set then partitions by
slot. Restricting the exact identity to slot `i` leaves only that part, which
needs at least the declared isolated rank `r_i`. Summing over the disjoint
parts gives `R >= sum_i r_i`; therefore positive fusion gain is impossible in
this construction class without mixed-product contact.
