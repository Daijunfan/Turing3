from __future__ import annotations


METHODS = (
    "independent evaluation", "Kronecker composition", "axis concatenation",
    "serendipitous substitution", "Pan pair fusion", "projection/peeling",
    "known direct-sum/aggregation", "best same-field explicit public certificate",
)


def baseline_envelope(candidates: list[dict]) -> dict:
    exact = [candidate for candidate in candidates if candidate.get("exact")]
    if not exact:
        return {"cost": None, "winner": None, "methods_required": list(METHODS)}
    winner = min(exact, key=lambda candidate: candidate["rank"])
    return {"cost": winner["rank"], "winner": winner, "methods_required": list(METHODS),
            "candidates": candidates}
