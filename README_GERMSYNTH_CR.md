# GermSynth-CR

Certificate-First Contact-Repair Algorithmogenesis audits inconsistent
external matrix-multiplication artifacts before allowing them into synthesis.

The main finding is forensic: FMM-Lille's current `<8,27,30>` page, Tensor,
LRP, and raw files represent ranks 3744, 3736, 3736, and 3825 respectively.
Tensor/LRP are identical and invalid; raw is a valid older rank-3825 program.
The recovered 21-term residual has flattening lower bound 9, disproving the
proposed eight-term truncation repair.

The trusted rank-250 base has six shared-U direction pairs but no shared-V
pairs, so the stated six `<4,3,3>:29` substitutions cannot reconstruct 3744
from the pinned inputs. The compiler retains the verified rank-3750
Kronecker fallback and quarantines all invalid candidates.

See:

- `results/latest_summary.json`
- `results/residual_autopsy.json`
- `results/raw_3825_explanation.json`
- `theory/contact_repair_model.md`
- `UPSTREAM_DISCREPANCY_REPORT_DRAFT.md`
