from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from typing import Any, Iterable

from .scheme_io import Scheme, SparseRow


def _mod_value(value: Fraction, modulus: int) -> int:
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def target_value(shape: tuple[int, int, int], a: int, b: int, c: int) -> int:
    m, n, p = shape
    ai, ak = divmod(a, n)
    bk, bj = divmod(b, p)
    ci, cj = divmod(c, p)
    return int(ai == ci and ak == bk and bj == cj)


def _residual_for_outputs(
    scheme: Scheme, outputs: Iterable[int], modulus: int | None = None, limit: int | None = None
) -> tuple[int, list[dict[str, Any]]]:
    wanted = set(outputs)
    accum: dict[int, dict[tuple[int, int], Fraction | int]] = {c: {} for c in wanted}
    for U, V, W in zip(scheme.U, scheme.V, scheme.W):
        for c, wc in W.items():
            if c not in wanted:
                continue
            out = accum[c]
            if modulus is None:
                for a, uc in U.items():
                    uw = uc * wc
                    for b, vc in V.items():
                        key = (a, b)
                        value = out.get(key, Fraction(0)) + uw * vc
                        if value:
                            out[key] = value
                        else:
                            out.pop(key, None)
            else:
                wv = _mod_value(wc, modulus)
                for a, uc in U.items():
                    uw = _mod_value(uc, modulus) * wv % modulus
                    for b, vc in V.items():
                        key = (a, b)
                        value = (int(out.get(key, 0)) + uw * _mod_value(vc, modulus)) % modulus
                        if value:
                            out[key] = value
                        else:
                            out.pop(key, None)
    residuals: list[dict[str, Any]] = []
    count = 0
    m, n, p = scheme.shape
    for c in sorted(wanted):
        actual = accum[c]
        ci, cj = divmod(c, p)
        for k in range(n):
            key = (ci * n + k, k * p + cj)
            expected = 1 if modulus is None else 1 % modulus
            value = actual.pop(key, 0)
            delta = value - expected if modulus is None else (int(value) - expected) % modulus
            if delta:
                count += 1
                if limit is None or len(residuals) < limit:
                    residuals.append({"a": key[0], "b": key[1], "c": c, "value": str(delta)})
        for (a, b), value in actual.items():
            if value:
                count += 1
                if limit is None or len(residuals) < limit:
                    residuals.append({"a": a, "b": b, "c": c, "value": str(value)})
    return count, residuals


def verify_scheme(scheme: Scheme, screen_primes: tuple[int, ...] = (1000003, 1000033)) -> dict[str, Any]:
    scheme.validate_dimensions()
    output_count = scheme.shape[0] * scheme.shape[2]
    field_modulus = 2 if scheme.field in ("GF(2)", "F2") else 3 if scheme.field in ("GF(3)", "F3") else None
    if field_modulus is not None:
        count, sample = _residual_for_outputs(scheme, range(output_count), field_modulus, 64)
        return {
            "status": "PASS" if count == 0 else "FAIL", "exact": count == 0,
            "exact_domain": f"GF({field_modulus})", "shape": list(scheme.shape), "rank": scheme.rank,
            "logical_tensor_coordinates": (scheme.shape[0] * scheme.shape[1])
                                          * (scheme.shape[1] * scheme.shape[2])
                                          * (scheme.shape[0] * scheme.shape[2]),
            "residual_count": count, "residual_sample": sample,
        }
    screens = []
    for prime in screen_primes:
        count, sample = _residual_for_outputs(scheme, range(output_count), prime, 8)
        screens.append({"prime": prime, "residual_count": count, "sample": sample})
        if count:
            return {"status": "FAIL", "exact": False, "modular_screens": screens,
                    "reason": f"nonzero residual modulo {prime}"}
    count, sample = _residual_for_outputs(scheme, range(output_count), None, 64)
    return {
        "status": "PASS" if count == 0 else "FAIL",
        "exact": count == 0,
        "shape": list(scheme.shape),
        "rank": scheme.rank,
        "coordinates_checked": scheme.shape[0] * scheme.shape[1] * scheme.shape[1] * scheme.shape[2]
                               * scheme.shape[0] * scheme.shape[2],
        "logical_tensor_coordinates": (scheme.shape[0] * scheme.shape[1])
                                      * (scheme.shape[1] * scheme.shape[2])
                                      * (scheme.shape[0] * scheme.shape[2]),
        "residual_count": count,
        "residual_sample": sample,
        "modular_screens": screens,
    }


def _normalize_product(U: SparseRow, V: SparseRow, W: SparseRow) -> tuple | None:
    if not U or not V or not W:
        return None
    def norm(row: SparseRow) -> tuple[tuple[tuple[int, Fraction], ...], Fraction]:
        pivot = row[min(row)]
        return tuple((index, value / pivot) for index, value in sorted(row.items())), pivot
    nu, su = norm(U)
    nv, sv = norm(V)
    nw, sw = norm(W)
    return nu, nv, nw, su * sv * sw


def audit_products(scheme: Scheme) -> dict[str, Any]:
    zero = []
    groups: dict[tuple, list[tuple[int, Fraction]]] = defaultdict(list)
    for product, factors in enumerate(zip(scheme.U, scheme.V, scheme.W)):
        normalized = _normalize_product(*factors)
        if normalized is None:
            zero.append(product)
            continue
        groups[normalized[:3]].append((product, normalized[3]))
    proportional = [members for members in groups.values() if len(members) > 1]
    mergeable = [members for members in proportional if sum(scale for _, scale in members) != 0]
    cancelling = [members for members in proportional if sum(scale for _, scale in members) == 0]
    return {
        "zero_product_ids": zero,
        "proportional_groups": [[[index, str(scale)] for index, scale in members] for members in proportional],
        "directly_mergeable_groups": [[[index, str(scale)] for index, scale in members] for members in mergeable],
        "cancelling_groups": [[[index, str(scale)] for index, scale in members] for members in cancelling],
    }


def tensor_hash(scheme: Scheme) -> str:
    payload = {"shape": scheme.shape, "field": scheme.field, "rank": scheme.rank,
               "U": [[(i, str(v)) for i, v in sorted(r.items())] for r in scheme.U],
               "V": [[(i, str(v)) for i, v in sorted(r.items())] for r in scheme.V],
               "W": [[(i, str(v)) for i, v in sorted(r.items())] for r in scheme.W]}
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()
