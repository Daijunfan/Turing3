# Residual completion

Let a target trilinear tensor be `T` and let a partial bilinear scheme be

`P = sum_i u_i tensor v_i tensor w_i`.

Its exact residual is `E(T,P) = T-P`. The residual completion rank is

`kappa(T,P) = rank(E(T,P))`.

This is an ordinary tensor-rank definition applied to a new workflow object;
it does not assert a new property of tensor rank.

## Machine representation

`Residual` stores every nonzero `(a,b,c,numerator,denominator)` coordinate in
canonical lexicographic JSONL. The header contains the shape, count, and a
SHA-256 over canonical coordinate lines. `compute_residual` expands all
rank-one products and checks the full target support, including all zeros.

For every factor axis `i`, flattening `E` into a matrix gives

`rank(E) >= rank(flatten_i(E))`.

The verifier computes all three ranks by exact sparse Gaussian elimination.
For the recovered 21-product component the ranks are `(6,9,6)`, so
`kappa >= 9`. Consequently an eight-term completion is impossible. This is a
general tensor-rank lower bound for this particular exact residual, not an
UNSAT claim from a restricted coefficient search.

## Local repair theorem

Let `A0`, `B0`, and `C0` be declared coordinate subspaces. Suppose the support
of `E(T,P)` is contained in `A0 × B0 × C0`, and let `D` be an exact rank-one
decomposition of `E` using factors supported in those subspaces. Then

`P + D = T`,

and every tensor coordinate outside `A0 × B0 × C0` is unchanged.

Proof: outside the product subspace, every rank-one term of `D` has at least
one zero factor coordinate, hence contributes zero. Inside it, `D=E=T-P` by
hypothesis. The two coordinate regions exhaust the tensor. The function
`local_repair_check` machine-checks the support premise, while
`verify_completion` checks `D=E` coordinate by coordinate.

## Replacement

Given a valid scheme `D` and removed term indices `S`, form

`P_S = D - sum_{i in S} q_i`.

Any exact completion of `T-P_S` may replace the removed terms. A strict rank
improvement occurs only when the completion rank is smaller than `|S|`.
Finite-field candidates must be lifted and reverified over the claimed field;
a finite-coefficient UNSAT result applies only to that finite search space.
