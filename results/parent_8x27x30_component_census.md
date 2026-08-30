# <8,27,30>:3736 component census

Overall: **FAIL** (`SOURCE_INVALID`)

The extractor recovered 238 ordinary slots and six two-slot candidates. Their product counts sum to 3736; their nominal savings sum to 14. Every exceptional component has nonzero exact residual, so the nominal 14 is not a proved fusion gain.

| Slots | Products | Nominal gain | Exact residuals | Standalone |
|---|---:|---:|---:|---|
| [19, 238] | 29 | 1 | 2232 | False |
| [23, 129] | 29 | 1 | 1728 | False |
| [25, 202] | 29 | 1 | 900 | False |
| [42, 48] | 29 | 1 | 21168 | False |
| [138, 175] | 29 | 1 | 1260 | False |
| [187, 223] | 21 | 9 | 936 | False |

Direct counterexample: `A[1,25] * B[25,10] -> C[1,10]` has expected coefficient 1 and actual coefficient 0.
