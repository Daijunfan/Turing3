from __future__ import annotations

from dataclasses import dataclass

from .exact_tensor import Scheme, verify_scheme
from .residual import Residual, compute_residual
from .residual_localization import flattening_ranks, localize_residual


@dataclass
class PartialScheme:
    target_shape: tuple[int, int, int]
    scheme: Scheme


def flattening_lower_bounds(residual: Residual):
    ranks = flattening_ranks(residual)
    return {"mode_ranks": list(ranks), "lower_bound": max(ranks, default=0)}


def verify_completion(residual: Residual, factors) -> bool:
    values = dict(residual.values)
    for U, V, W in factors:
        for a, uc in U.items():
            for b, vc in V.items():
                for c, wc in W.items():
                    key = (a, b, c)
                    value = values.get(key, 0) - uc * vc * wc
                    if value:
                        values[key] = value
                    else:
                        values.pop(key, None)
    return not values


def replace_terms(scheme: Scheme, removed_terms: set[int], completion) -> Scheme:
    kept = [index for index in range(scheme.rank) if index not in removed_terms]
    products = [(scheme.U[index], scheme.V[index], scheme.W[index]) for index in kept] + list(completion)
    return Scheme(scheme.field, scheme.shape, len(products), [p[0] for p in products], [p[1] for p in products],
                  [p[2] for p in products], "residual replacement", {"removed_terms": sorted(removed_terms)})


def local_repair_check(residual: Residual, allowed_support: tuple[set[int], set[int], set[int]]) -> dict:
    outside = [coordinate for coordinate in residual.values
               if any(coordinate[axis] not in allowed_support[axis] for axis in range(3))]
    return {"pass": not outside, "outside_support_count": len(outside), "outside_support_sample": outside[:32],
            "theorem": "A completion supported inside the declared factor subspaces changes no outside tensor coordinate."}


def search_completion(residual: Residual, rank_budget: int) -> dict:
    lower = max(flattening_ranks(residual), default=0)
    if lower > rank_budget:
        return {"status": "PROVED_IMPOSSIBLE_BY_FLATTENING", "rank_budget": rank_budget, "lower_bound": lower}
    return {"status": "UNKNOWN", "rank_budget": rank_budget, "lower_bound": lower,
            "reason": "No unrestricted tensor-rank solver completed."}


def emit_certificate(residual: Residual, completion, source: str) -> dict:
    return {"format": "germsynth-cr-residual-completion-v1", "source": source,
            "residual_sha256": residual.sha256(), "rank": len(completion),
            "verified": verify_completion(residual, completion)}
