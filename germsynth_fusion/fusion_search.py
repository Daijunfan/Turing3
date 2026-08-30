from __future__ import annotations

import json
from pathlib import Path
from typing import Any


PRIORITY_SHAPES = [(2, 2, 2), (2, 2, 3), (2, 3, 3), (2, 4, 4), (3, 3, 3)]


def write_checkpoint(path: str | Path, result: dict[str, Any]) -> None:
    Path(path).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


def source_gated_search(parent_exact: bool, extracted_kernels: list[dict[str, Any]], path: str | Path) -> dict[str, Any]:
    """Record a deterministic search gate: invalid seeds are never promoted as discoveries."""
    result = {
        "status": "NOT_RUN" if not parent_exact else "COMPLETED",
        "reason": "SOURCE_INVALID: exact parent tensor verification failed" if not parent_exact else None,
        "priority_shapes": [list(shape) for shape in PRIORITY_SHAPES],
        "valid_seed_count": sum(kernel.get("status") == "PASS" for kernel in extracted_kernels),
        "novel_fusion_kernel": False,
        "new_exact_rank": False,
        "search_representation": "exact sparse rank-one terms with checkpointing",
        "attempts": [],
    }
    write_checkpoint(path, result)
    return result
