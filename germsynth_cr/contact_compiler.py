from __future__ import annotations

from collections import defaultdict
from fractions import Fraction

from .contact_hypergraph import ContactHypergraph, ContactPackingSolver
from .contact_rank import ContactClass, ContactOrgan
from .exact_tensor import Scheme, verify_scheme
from germsynth_fusion.fusion_constructor import kronecker_scheme


def _direction(row):
    pivot = row[min(row)]
    return tuple((index, value / pivot) for index, value in sorted(row.items()))


def proportional_pairs(scheme: Scheme, factor: str) -> list[tuple[int, int]]:
    rows = getattr(scheme, factor)
    groups = defaultdict(list)
    for product, row in enumerate(rows):
        groups[_direction(row)].append(product)
    return [tuple(group) for group in groups.values() if len(group) == 2]


def analyze_rank3744_preconditions(outer: Scheme, inner: Scheme, concat: Scheme) -> dict:
    groups = {factor: proportional_pairs(outer, factor) for factor in ("U", "V", "W")}
    # A <4,3,3> concatenation of two <2,3,3> slots needs a shared V direction:
    # local factor dimensions become 12,9,12. Shared U instead gives 6,18,12 (<2,3,6>).
    required = "V"
    available = groups[required]
    return {
        "factor_direction_pairs": {factor: [list(pair) for pair in pairs] for factor, pairs in groups.items()},
        "required_shared_factor_for_4x3x3": required,
        "available_pair_count": len(available),
        "preconditions_met": len(available) >= 6 and concat.shape == (4, 3, 3) and concat.rank == 29,
        "shared_U_local_factor_dimensions": [6, 18, 12],
        "required_4x3x3_factor_dimensions": [12, 9, 12],
    }


class ContactCompiler:
    def compile_rank3744(self, outer: Scheme, inner: Scheme, concat: Scheme):
        analysis = analyze_rank3744_preconditions(outer, inner, concat)
        if not analysis["preconditions_met"]:
            return {"status": "FAIL", "reason": "trusted outer scheme has no six shared-V pairs",
                    "analysis": analysis, "fallback": kronecker_scheme(outer, inner)}
        raise NotImplementedError("shared-V embedding is not exercised by the pinned outer certificate")

    def verify(self, scheme: Scheme):
        return verify_scheme(scheme)
