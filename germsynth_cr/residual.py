from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from .exact_tensor import Scheme


Coordinate = tuple[int, int, int]


@dataclass
class Residual:
    shape: tuple[int, int, int]
    values: dict[Coordinate, Fraction]

    def sha256(self) -> str:
        digest = hashlib.sha256()
        for coordinate, value in sorted(self.values.items()):
            digest.update(f"{coordinate[0]},{coordinate[1]},{coordinate[2]}:{value}\n".encode())
        return digest.hexdigest()

    def write_jsonl(self, path: str | Path) -> None:
        with Path(path).open("w") as output:
            header = {"format": "germsynth-cr-residual-jsonl-v1", "shape": self.shape,
                      "nonzero_count": len(self.values), "sha256": self.sha256()}
            output.write(json.dumps(header, separators=(",", ":")) + "\n")
            for (a, b, c), value in sorted(self.values.items()):
                output.write(json.dumps([a, b, c, value.numerator, value.denominator], separators=(",", ":")) + "\n")


def compute_residual(target_shape: tuple[int, int, int], partial: Scheme) -> Residual:
    """Return E=T_mm-P exactly, retaining every and only nonzero coordinate."""
    if partial.shape != target_shape:
        raise ValueError("target and partial shapes differ")
    m, n, p = target_shape
    by_output: dict[int, dict[tuple[int, int], Fraction]] = {c: {} for c in range(m * p)}
    for U, V, W in zip(partial.U, partial.V, partial.W):
        for c, wc in W.items():
            output = by_output[c]
            for a, uc in U.items():
                for b, vc in V.items():
                    key = (a, b)
                    value = output.get(key, Fraction()) - uc * vc * wc
                    if value:
                        output[key] = value
                    else:
                        output.pop(key, None)
    for i in range(m):
        for k in range(n):
            for j in range(p):
                c, key = i * p + j, (i * n + k, k * p + j)
                value = by_output[c].get(key, Fraction()) + 1
                if value:
                    by_output[c][key] = value
                else:
                    by_output[c].pop(key, None)
    return Residual(target_shape, {(a, b, c): value for c, output in by_output.items()
                                   for (a, b), value in output.items() if value})


def residual_from_factors(shape: tuple[int, int, int], products) -> Residual:
    scheme = Scheme("Q", shape, len(products), [p[0] for p in products], [p[1] for p in products],
                    [p[2] for p in products], "partial")
    return compute_residual(shape, scheme)
