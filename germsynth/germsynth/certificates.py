"""Serialization helpers for independently checkable germ certificates."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .germ import BilinearGerm


def write_certificate(
    germ: BilinearGerm,
    path: str | Path,
    *,
    experiments: dict[str, Any] | None = None,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = germ.certificate_payload()
    payload["proof_obligation"] = {
        "local": "target tensor equals the exact sum of rank-one term tensors",
        "recursive_closure": (
            "For each power-of-two scale, substitute blocks for scalar variables; "
            "the local polynomial identity and the induction hypothesis imply the next scale."
        ),
    }
    if experiments is not None:
        payload["experiments"] = experiments
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return destination


def read_certificate(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))
