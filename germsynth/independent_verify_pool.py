#!/usr/bin/env python3
"""Independent standard-library verifier for the polymorphic germ pool."""
from __future__ import annotations

from itertools import combinations
import json
from pathlib import Path
import sys


def target_tensor() -> tuple[int, ...]:
    target = [0] * 64
    for i in range(2):
        for j in range(2):
            for k in range(2):
                a, b, c = i * 2 + k, k * 2 + j, i * 2 + j
                target[(a * 4 + b) * 4 + c] += 1
    return tuple(target)


def exact_tensor(terms: list[dict[str, list[int]]]) -> tuple[int, ...]:
    result = [0] * 64
    for term in terms:
        u, v, w = term["u"], term["v"], term["w"]
        if len(u) != 4 or len(v) != 4 or len(w) != 4:
            raise AssertionError("invalid term dimensions")
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    result[(i * 4 + j) * 4 + k] += u[i] * v[j] * w[k]
    return tuple(result)


def support_mask(values: list[int]) -> int:
    return sum((1 << index) for index, value in enumerate(values) if value)


def resources(terms: list[dict[str, list[int]]]) -> frozenset[tuple[str, int]]:
    result = set()
    for term in terms:
        result.add(("left", support_mask(term["u"])))
        result.add(("right", support_mask(term["v"])))
        result.add(("output", support_mask(term["w"])))
    return frozenset(result)


def coverage(requirements, universe, order, selected):
    total = covered = 0
    for failed in combinations(universe, order):
        total += 1
        if any(requirements[index].isdisjoint(failed) for index in selected):
            covered += 1
    return covered, total


def main() -> None:
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "certificates/matrix_phenotype_pool.json")
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != "germsynth-phenotype-pool-v1":
        raise AssertionError("unsupported pool schema")
    phenotypes = data["phenotypes"]
    if len(phenotypes) != data["phenotype_count"]:
        raise AssertionError("phenotype count mismatch")
    target = target_tensor()
    requirements = []
    coefficient_keys = set()
    for phenotype in phenotypes:
        terms = phenotype["terms"]
        if len(terms) != data["rank"]:
            raise AssertionError("phenotype rank mismatch")
        if exact_tensor(terms) != target:
            raise AssertionError(f"invalid phenotype: {phenotype['name']}")
        key = tuple(tuple(term[mode]) for term in terms for mode in ("u", "v", "w"))
        coefficient_keys.add(key)
        requirements.append(resources(terms))
    if len(coefficient_keys) != len(phenotypes):
        raise AssertionError("pool contains duplicate coefficient phenotypes")
    universe = tuple(sorted(set().union(*requirements)))

    for order in (1, 2, 3):
        actual = coverage(requirements, universe, order, tuple(range(len(phenotypes))))
        stated = data["failure_coverage"][str(order)]
        if actual != (stated["pool_survives"], stated["failure_sets"]):
            raise AssertionError(f"coverage mismatch for order {order}")

    cover_data = data["minimum_single_fault_cover"]
    selected = tuple(cover_data["phenotype_indices"])
    selected_coverage = coverage(requirements, universe, 1, selected)
    if selected_coverage[0] != selected_coverage[1]:
        raise AssertionError("stated regenerative basis does not cover all single failures")
    # Prove minimality by exhaustive rejection of every smaller subset.
    for size in range(1, len(selected)):
        for candidate in combinations(range(len(phenotypes)), size):
            covered, total = coverage(requirements, universe, 1, candidate)
            if covered == total:
                raise AssertionError("a smaller full single-fault cover exists")
    if coverage(requirements, universe, 2, selected) != (
        cover_data["pairs"]["covered"], cover_data["pairs"]["total"]
    ):
        raise AssertionError("selected-pool pair coverage mismatch")
    if coverage(requirements, universe, 3, selected) != (
        cover_data["triples"]["covered"], cover_data["triples"]["total"]
    ):
        raise AssertionError("selected-pool triple coverage mismatch")

    print(json.dumps({
        "status": "PASS",
        "certificate": str(path),
        "exact_phenotypes": len(phenotypes),
        "resource_universe": len(universe),
        "minimum_single_fault_cover_size": len(selected),
        "single_failures_covered": selected_coverage[0],
        "single_failures_total": selected_coverage[1],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
