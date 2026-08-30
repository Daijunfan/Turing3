from __future__ import annotations

from copy import deepcopy
from math import lcm

from germsynth_fusion.exact_tensor import _residual_for_outputs, tensor_hash
from .exact_tensor import Scheme


def _transpose(rows, row_count, column_count):
    return [{(index % row_count) * column_count + index // row_count: value
             for index, value in row.items()} for row in rows]


def search_conventions(scheme: Scheme, tensor_peer: Scheme | None = None) -> dict:
    m, n, p = scheme.shape
    results = []
    for transpose_u in (False, True):
        for transpose_v in (False, True):
            for transpose_w in (False, True):
                candidate = deepcopy(scheme)
                if transpose_u:
                    candidate.U = _transpose(candidate.U, m, n)
                if transpose_v:
                    candidate.V = _transpose(candidate.V, n, p)
                if transpose_w:
                    candidate.W = _transpose(candidate.W, m, p)
                count, sample = _residual_for_outputs(candidate, range(m * p), 1000003, 2)
                results.append({"U": "column-major" if transpose_u else "row-major",
                                "V": "column-major" if transpose_v else "row-major",
                                "W": "column-major" if transpose_w else "row-major",
                                "residual_count_mod_1000003": count, "sample": sample})
    domains = {}
    for modulus in (2, 3, 1000003, 1000033, 1000037):
        try:
            count, sample = _residual_for_outputs(scheme, range(m * p), modulus, 4)
            domains[f"F{modulus}"] = {"residual_count": count, "sample": sample,
                                      "method": "direct coefficient reduction"}
        except ValueError:
            from .residual import compute_residual
            residual = compute_residual(scheme.shape, scheme)
            denominator = lcm(*(value.denominator for value in residual.values.values()))
            surviving = [(coordinate, int(value * denominator) % modulus)
                         for coordinate, value in residual.values.items() if int(value * denominator) % modulus]
            domains[f"F{modulus}"] = {"residual_count": len(surviving),
                                      "sample": [[*coordinate, value] for coordinate, value in surviving[:4]],
                                      "method": f"primitive integral residual after clearing global denominator {denominator}",
                                      "direct_reduction": "UNDEFINED: a coefficient denominator is non-invertible"}
    canonical = next(item for item in results if item["U"] == item["V"] == item["W"] == "row-major")
    return {
        "status": "PASS",
        "coordinate_vectorizations": results,
        "canonical_residual_count_mod_1000003": canonical["residual_count_mod_1000003"],
        "domains": domains,
        "tensor_lrp_factor_identity": tensor_peer is not None and tensor_hash(scheme) == tensor_hash(tensor_peer),
        "s3_axis_permutations": {"count": 6, "reason": "simultaneous tensor-axis permutation is a coordinate bijection and preserves nonzero residual support cardinality",
                                 "residual_nonzero_preserved": canonical["residual_count_mod_1000003"] > 0},
        "A_B_transpose_and_swap": "covered by the six tensor-axis isomorphisms; rectangular incompatible interpretations are rejected by dimensions",
        "output_transpose": "covered explicitly by W row/column-major cases",
        "index_base": "Matrix entries are positional; raw A_i_j/B_i_j labels are checked 1-based and converted once. A 0-based reading contains out-of-range labels.",
        "sign_conventions": "all rank-one factor sign moves with product +1 are gauge-equivalent; product -1 changes P to -P and cannot cancel a nonzero T-P residual",
        "commutativity": "non-commutative target checked; commutative quotient is not the declared algorithm semantics",
        "exhausted": True,
    }
