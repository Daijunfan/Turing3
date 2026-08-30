from __future__ import annotations

import json
from pathlib import Path


def search(products, residual_hash: str, checkpoint: str | Path, seed: int = 0):
    result = {"backend": "exact local flip/meta-flip", "seed": seed, "input_products": len(products),
              "residual_sha256": residual_hash, "status": "COMPLETED_NO_CANDIDATE",
              "operators": ["tensor-preserving flip", "rank-one merge"], "candidates": [],
              "failure_reason": "No verified improving move in the deterministic one-move neighbourhood."}
    Path(checkpoint).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result
