# HANDOFF TO BRAIN — GermSynth-F

## Status table

| Gate | Status |
|---|---|
| BASELINE_REPRODUCTION | **PASS** |
| PARENT_EXACT_VERIFICATION | **FAIL** |
| AUTOMATIC_SLOT_INFERENCE | **PASS** |
| FUSION_KERNEL_EXTRACTION | **FAIL** |
| FULL_GAIN_14_EXPLAINED | **FAIL** |
| CXX_INDEPENDENT_VERIFICATION | **FAIL** |
| CONSTRUCTOR_PASS | **FAIL** |
| ALL_SIX_DEFICITS_CLOSED | **FAIL** |
| NOVEL_FUSION_KERNEL | **FAIL** |
| NEW_EXACT_RANK | **FAIL** |
| NEW_EXPONENT | **FAIL** |
| TURING_PATH_PASS | **FAIL** |

## Three most important findings

1. FMM-Lille is internally inconsistent: Triad/LRP declare 3736, the current page declares 3744, and the raw evaluation file contains 3825 products.
2. Exact verification rejects that artifact: 28,098 tensor coordinates are nonzero modulo 1,000,003; a direct rational counterexample is `A[1,25] B[25,10] -> C[1,10]`, expected 1, actual 0.
3. Automatic extraction recovers 238 ordinary slots plus six two-slot candidates (five rank 29 and one rank 21), whose nominal savings are 5×1+9=14, but every candidate has a nonzero exact residual.

## Known versus new

The 29-for-30 narrative and the stale 3736 discussion are pre-existing claims in the downloaded catalog draft. GermSynth-F's new result here is a reproducible falsification of the supplied artifact, not a new fast multiplication identity.

## The nominal 14-product explanation

`238×15 + 5×29 + 1×21 = 3736`, versus `250×15 = 3750`. This is only rank arithmetic. Exact residuals disprove it as a certificate for the downloaded artifact.

## Is the rank-29 kernel standalone?

No. All five rank-29 candidates have nonzero residuals; the representative certificate is deliberately marked FAIL.

## Constructor universality

Not established. The constructor rejects every extracted candidate because none has a PASS certificate.

## New rank or exponent

None. Single-type spectral regressions recover log₂3 and log₂7, but invalid source data prevents a fusion before/after exponent comparison.

## Failed hypotheses and counterexamples

- Hypothesis: the downloadable `<8,27,30>:3736` Triads exactly compute matrix multiplication. Counterexample above.
- Hypothesis: a recovered 29-product pair is standalone. Each of five candidates has nonzero residual (see component census).
- Hypothesis: rank arithmetic alone explains a valid 14 gain. It does not; the full exceptional closure retains the parent's 28,098 modular residuals.

## Three next research questions

1. Can FMM-Lille supply a corrected explicit artifact matching either rank 3744 or the historical 3736 claim?
2. Which generation step introduced the six invalid exceptional components, and can the missing exact terms be reconstructed from provenance?
3. Do any verified external fusion artifacts yield a standalone two-slot kernel under the same exact extractor?

## One-command reproduction

```bash
./reproduce_fusion.sh
```

## Formal no-contact proposition

Within non-overlap composition, if every product is incident to at most one slot, product sets partition by slot. Exact restriction to each slot therefore needs at least that slot's declared rank; summing gives R >= sum_i r_i, so fusion_gain <= 0.
