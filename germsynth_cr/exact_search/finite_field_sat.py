from __future__ import annotations

import json
import time
from pathlib import Path


def search(residual, rank_budget: int, field: int, checkpoint: str | Path, seed: int = 0):
    started = time.time()
    result = {"backend": "finite-field SAT", "field": field, "rank_budget": rank_budget, "seed": seed,
              "search_space": "Brent equations over the declared residual support",
              "status": "NOT_RUN_NO_SAT_SOLVER", "scope": "No UNSAT claim", "seconds": time.time() - started,
              "candidates": []}
    Path(checkpoint).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result
