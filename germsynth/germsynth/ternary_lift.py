"""Exact support-first lift from GF(2) masks to ternary integer coefficients.

For a fixed GF(2) support, signs in each rank-one term are the only unknowns.
Two sign gauges are removed per term.  A meet-in-the-middle exact tensor sum
then decides whether the support admits a {-1,0,1} realization over integers.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import prod
from time import perf_counter
from typing import Iterable

from .gf2_flip import MaskScheme
from .tensor import RankOneTerm, matrix_multiplication_2x2_specification


@dataclass(frozen=True, slots=True)
class LiftResult:
    terms: tuple[RankOneTerm, ...]
    partition_left: tuple[int, ...]
    partition_right: tuple[int, ...]
    left_assignments: int
    right_assignments: int
    unique_left_sums: int
    elapsed_seconds: float


def _support(mask: int, dimension: int = 4) -> tuple[int, ...]:
    return tuple(index for index in range(dimension) if (mask >> index) & 1)


def signed_variants(mask_term: tuple[int, int, int]) -> tuple[RankOneTerm, ...]:
    supports = tuple(_support(mask) for mask in mask_term)
    if any(not support for support in supports):
        raise ValueError("rank-one support must be nonempty in every mode")
    u_support, v_support, w_support = supports
    free = (
        tuple(("u", index) for index in u_support[1:])
        + tuple(("v", index) for index in v_support[1:])
        + tuple(("w", index) for index in w_support)
    )
    variants: list[RankOneTerm] = []
    for bits in range(1 << len(free)):
        u, v, w = [0] * 4, [0] * 4, [0] * 4
        # Two tensor sign gauges are fixed by making the first nonzero u/v entries +1.
        u[u_support[0]] = 1
        v[v_support[0]] = 1
        for bit, (mode, index) in enumerate(free):
            value = -1 if (bits >> bit) & 1 else 1
            if mode == "u":
                u[index] = value
            elif mode == "v":
                v[index] = value
            else:
                w[index] = value
        variants.append(RankOneTerm(tuple(u), tuple(v), tuple(w)))
    return tuple(variants)


def _balanced_partition(variant_counts: list[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Choose a partition minimizing the larger Cartesian-product side."""
    indices = tuple(range(len(variant_counts)))
    best: tuple[int, int, tuple[int, ...], tuple[int, ...]] | None = None
    # Pin index 0 to the left to remove symmetric duplicates.
    for left_size in range(1, len(indices)):
        for left_tail in combinations(indices[1:], left_size - 1):
            left = (0,) + left_tail
            left_set = set(left)
            right = tuple(index for index in indices if index not in left_set)
            left_product = prod(variant_counts[index] for index in left)
            right_product = prod(variant_counts[index] for index in right)
            key = (max(left_product, right_product), left_product + right_product, left, right)
            if best is None or key < best:
                best = key
    if best is None:
        raise ValueError("at least two terms are required for meet-in-the-middle")
    return best[2], best[3]


def _enumerate_sums(
    variant_lists: list[tuple[RankOneTerm, ...]], indices: tuple[int, ...], tensor_size: int
) -> tuple[dict[tuple[int, ...], tuple[RankOneTerm, ...]], int]:
    sums: dict[tuple[int, ...], tuple[RankOneTerm, ...]] = {}
    generated = 0

    def visit(position: int, current: list[int], choice: list[RankOneTerm]) -> None:
        nonlocal generated
        if position == len(indices):
            generated += 1
            sums.setdefault(tuple(current), tuple(choice))
            return
        term_index = indices[position]
        for term in variant_lists[term_index]:
            tensor = term.tensor()
            for coordinate, value in enumerate(tensor):
                current[coordinate] += value
            choice.append(term)
            visit(position + 1, current, choice)
            choice.pop()
            for coordinate, value in enumerate(tensor):
                current[coordinate] -= value

    visit(0, [0] * tensor_size, [])
    return sums, generated


def lift_support_to_integers(scheme: MaskScheme) -> LiftResult:
    start_time = perf_counter()
    specification = matrix_multiplication_2x2_specification()
    variant_lists = [signed_variants(term) for term in scheme]
    counts = [len(variants) for variants in variant_lists]
    left_indices, right_indices = _balanced_partition(counts)
    left_sums, left_generated = _enumerate_sums(
        variant_lists, left_indices, len(specification.tensor)
    )

    answer: tuple[RankOneTerm, ...] | None = None
    right_generated = 0

    def visit_right(position: int, current: list[int], choice: list[RankOneTerm]) -> None:
        nonlocal answer, right_generated
        if answer is not None:
            return
        if position == len(right_indices):
            right_generated += 1
            needed = tuple(
                target - partial
                for target, partial in zip(specification.tensor, current, strict=True)
            )
            left_choice = left_sums.get(needed)
            if left_choice is not None:
                by_index: dict[int, RankOneTerm] = {}
                for index, term in zip(left_indices, left_choice, strict=True):
                    by_index[index] = term
                for index, term in zip(right_indices, choice, strict=True):
                    by_index[index] = term
                answer = tuple(by_index[index] for index in range(len(scheme)))
            return

        term_index = right_indices[position]
        for term in variant_lists[term_index]:
            tensor = term.tensor()
            for coordinate, value in enumerate(tensor):
                current[coordinate] += value
            choice.append(term)
            visit_right(position + 1, current, choice)
            choice.pop()
            for coordinate, value in enumerate(tensor):
                current[coordinate] -= value
            if answer is not None:
                return

    visit_right(0, [0] * len(specification.tensor), [])
    if answer is None:
        raise RuntimeError("GF(2) support has no ternary integer lift under the chosen gauges")

    reconstructed = [0] * len(specification.tensor)
    for term in answer:
        for index, value in enumerate(term.tensor()):
            reconstructed[index] += value
    if tuple(reconstructed) != specification.tensor:
        raise AssertionError("meet-in-the-middle lift failed exact tensor verification")

    return LiftResult(
        terms=answer,
        partition_left=left_indices,
        partition_right=right_indices,
        left_assignments=left_generated,
        right_assignments=right_generated,
        unique_left_sums=len(left_sums),
        elapsed_seconds=perf_counter() - start_time,
    )
