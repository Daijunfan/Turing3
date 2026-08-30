from __future__ import annotations

from ..residual_completion import search_completion


def solve(residual, rank_budget: int):
    return {"backend": "exact residual solver", **search_completion(residual, rank_budget)}
