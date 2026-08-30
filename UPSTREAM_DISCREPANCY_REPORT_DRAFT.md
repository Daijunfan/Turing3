# Draft upstream discrepancy report: FMM-Lille `<8,27,30>`

Observed on 2026-08-31:

| Representation | Declared/semantic rank | Exact result | Last-Modified |
|---|---:|---|---|
| Web page | 3744 | metadata only | 2026-06-01 21:52:15 GMT |
| Tensor TriadSet | 3736 | FAIL, 28,098 residual coordinates | 2026-06-01 21:52:16 GMT |
| LRP | 3736 | factor-identical to Tensor; same FAIL | 2026-06-01 21:52:16 GMT |
| raw evaluation | 3825 | PASS over Q | 2022-11-23 14:32:31 GMT |

The exact counterexample `A[1,25] B[25,10] -> C[1,10]` has expected
coefficient 1 and actual coefficient 0 in Tensor/LRP.

The rank-3736 artifact decomposes arithmetically as 238 ordinary 15-product
slots, five 29-product candidate components, and one 21-product component.
None of the six exceptional components is exact. The 21-product residual has
flattening ranks `(6,9,6)`, proving it cannot be repaired with eight rank-one
terms. Thus it is not simply a rank-3744 file truncated by eight terms.

The pinned rank-250 base has six proportional U-factor pairs and no
proportional V-factor pairs. Pairing two `<2,3,3>` slots through shared U gives
factor-space dimensions `(6,18,12)`, whereas `<4,3,3>` has `(12,9,12)`.
Consequently the page recipe cannot be reconstructed from the currently
published base certificate by the stated six `<4,3,3>:29` substitutions.

Please provide the exact rank-3744 factors and the precise rank-250 base
representative/orientation if a different valid artifact supports the page.
All hashes and machine-readable evidence are in `sources.lock.json` and
`results/residual_autopsy.json`.
