# GermSynth exact experiment results

**Overall status:** PASS

## Discovered germs

| Problem | Automatically discovered rank | Exact base cases | Recursive exponent |
|---|---:|---:|---:|
| Degree-1 polynomial multiplication | 3 | 625 | 1.584962500721 |
| 2x2 matrix multiplication | 7 | 6561 | 2.807354922058 |

The matrix search ran 20 deterministic seeds, found a rank-7 support and an exact integer lift in every run, and produced 16 distinct exact coefficient phenotypes.  The selected primary support was reached from schoolbook rank 8 after 171 exact local flips.

## Largest recursive checks

| Algorithm family | Largest size | Exact scalar multiplications | Schoolbook count | Equality |
|---|---:|---:|---:|---|
| Matrix multiplication | 128 x 128 | 823543 | 2097152 | PASS |
| Polynomial multiplication | 1024 coefficients | 59049 | 1048576 | PASS |

A 64x64 matrix product using a different exact phenotype at arbitrary recursive nodes also passed, with exactly 117649 scalar multiplications.  Under the stated abstract resource model, exact subset optimization found a minimal three-phenotype basis that survives all 44 single-resource failures; a strict 32x32 run sampled one such failure independently at every recursive node, allowed no fallback, and passed.  The full 16-phenotype pool survives 871 of 946 double-resource failures.

## Verification boundary

The local identities use exact integer arithmetic.  Universal power-of-two correctness follows by structural induction on block substitution.  Random large-size checks verify the implementation/proof correspondence; they are not the logical basis of the universal theorem.
