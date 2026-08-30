#!/usr/bin/env python3
"""Independent, standard-library verifier for GermSynth JSON certificates.

This file deliberately does not import the GermSynth package.  It reconstructs
problem tensors, rank-one tensors, hashes, exhaustive base evaluations and the
closed-form scalar multiplication count directly from JSON.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from math import isclose, log
from pathlib import Path
from typing import Sequence


def tensor_hash(values: Sequence[int]) -> str:
    return sha256(",".join(str(value) for value in values).encode("ascii")).hexdigest()


def matrix_target() -> tuple[int, ...]:
    target = [0] * 64
    for i in range(2):
        for j in range(2):
            for k in range(2):
                a = i * 2 + k
                b = k * 2 + j
                c = i * 2 + j
                target[(a * 4 + b) * 4 + c] += 1
    return tuple(target)


def polynomial_target() -> tuple[int, ...]:
    target = [0] * 12
    for i, j, k in ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 2)):
        target[(i * 2 + j) * 3 + k] += 1
    return tuple(target)


def outer(u: Sequence[int], v: Sequence[int], w: Sequence[int]) -> tuple[int, ...]:
    return tuple(a * b * c for a in u for b in v for c in w)


def evaluate_terms(terms: list[dict[str, list[int]]], a: Sequence[int], b: Sequence[int]) -> tuple[int, ...]:
    output_dimension = len(terms[0]["w"])
    output = [0] * output_dimension
    for term in terms:
        left = sum(coefficient * value for coefficient, value in zip(term["u"], a, strict=True))
        right = sum(coefficient * value for coefficient, value in zip(term["v"], b, strict=True))
        value = left * right
        for index, coefficient in enumerate(term["w"]):
            output[index] += coefficient * value
    return tuple(output)


def evaluate_target(target: Sequence[int], a_dim: int, b_dim: int, out_dim: int, a: Sequence[int], b: Sequence[int]) -> tuple[int, ...]:
    output = [0] * out_dim
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            for k in range(out_dim):
                output[k] += target[(i * b_dim + j) * out_dim + k] * ai * bj
    return tuple(output)


def verify(path: Path, exhaustive: bool = True) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != "germsynth-bilinear-germ-v1":
        raise AssertionError("unsupported certificate schema")
    problem = data["problem"]
    kind = problem["kind"]
    expected_target = matrix_target() if kind == "matrix_multiplication" else polynomial_target()
    target = tuple(data["target_tensor"])
    if target != expected_target:
        raise AssertionError("certificate target does not match the declared problem")
    if tensor_hash(target) != data["target_tensor_sha256"]:
        raise AssertionError("target tensor hash mismatch")

    a_dim, b_dim, out_dim = problem["a_dim"], problem["b_dim"], problem["out_dim"]
    terms = data["terms"]
    if len(terms) != data["rank"]:
        raise AssertionError("rank field does not match term count")
    reconstructed = [0] * len(target)
    for term in terms:
        if len(term["u"]) != a_dim or len(term["v"]) != b_dim or len(term["w"]) != out_dim:
            raise AssertionError("rank-one term has invalid dimensions")
        if any(value not in (-1, 0, 1) for mode in ("u", "v", "w") for value in term[mode]):
            raise AssertionError("prototype certificate is expected to be ternary")
        contribution = outer(term["u"], term["v"], term["w"])
        for index, value in enumerate(contribution):
            reconstructed[index] += value
    reconstructed_tuple = tuple(reconstructed)
    if reconstructed_tuple != target:
        raise AssertionError("exact tensor identity failed")
    if tensor_hash(reconstructed_tuple) != data["reconstructed_tensor_sha256"]:
        raise AssertionError("reconstructed tensor hash mismatch")
    if not data["local_identity_verified"]:
        raise AssertionError("certificate does not assert local verification")

    block_factor = data["block_factor"]
    expected_exponent = log(len(terms), block_factor)
    stated_exponent = data["complexity_certificate"]["scalar_multiplication_exponent"]
    if not isclose(expected_exponent, stated_exponent, rel_tol=0.0, abs_tol=1e-12):
        raise AssertionError("complexity exponent mismatch")

    exhaustive_cases = 0
    if exhaustive:
        alphabet = (-1, 0, 1) if kind == "matrix_multiplication" else (-2, -1, 0, 1, 2)
        for a in product(alphabet, repeat=a_dim):
            for b in product(alphabet, repeat=b_dim):
                exhaustive_cases += 1
                actual = evaluate_terms(terms, a, b)
                expected = evaluate_target(target, a_dim, b_dim, out_dim, a, b)
                if actual != expected:
                    raise AssertionError(f"exhaustive base evaluation failed at a={a}, b={b}")

    return {
        "certificate": str(path),
        "name": data["name"],
        "rank": len(terms),
        "block_factor": block_factor,
        "exponent": expected_exponent,
        "exact_tensor_identity": True,
        "exhaustive_base_cases": exhaustive_cases,
        "status": "PASS",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificates", nargs="+", type=Path)
    parser.add_argument("--skip-exhaustive", action="store_true")
    arguments = parser.parse_args()
    for certificate in arguments.certificates:
        result = verify(certificate, exhaustive=not arguments.skip_exhaustive)
        print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
