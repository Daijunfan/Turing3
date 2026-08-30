"""Exact local search in the GF(2) flip graph for 2x2 matrix multiplication.

The local flip/reduction identities are verified after every generated move.
The implementation starts from the rank-8 schoolbook tensor decomposition and
searches for a rank-7 support topology.  Coefficients are lifted to the integer
ring separately by :mod:`germsynth.ternary_lift`.
"""
from __future__ import annotations

from dataclasses import dataclass
import random
from time import perf_counter
from typing import Iterable

MaskTerm = tuple[int, int, int]
MaskScheme = tuple[MaskTerm, ...]

PERMUTATIONS: tuple[tuple[int, int, int], ...] = (
    (0, 1, 2),
    (0, 2, 1),
    (1, 0, 2),
    (1, 2, 0),
    (2, 0, 1),
    (2, 1, 0),
)


@dataclass(frozen=True, slots=True)
class FlipSearchResult:
    seed: int
    scheme: MaskScheme
    flips: int
    visited_states: int
    elapsed_seconds: float
    naive_xor_additions: int


def canonical_scheme(terms: Iterable[MaskTerm]) -> MaskScheme:
    return tuple(sorted(terms))


def tensor_of_term(term: MaskTerm) -> int:
    """Return the 4x4x4 GF(2) tensor as one 64-bit integer."""
    a, b, c = term
    result = 0
    for i in range(4):
        if (a >> i) & 1:
            for j in range(4):
                if (b >> j) & 1:
                    for k in range(4):
                        if (c >> k) & 1:
                            result ^= 1 << ((i * 4 + j) * 4 + k)
    return result


def tensor_of_scheme(scheme: MaskScheme) -> int:
    result = 0
    for term in scheme:
        result ^= tensor_of_term(term)
    return result


def target_matrix_multiplication_tensor() -> int:
    result = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ai = i * 2 + k
                bj = k * 2 + j
                ci = i * 2 + j
                result ^= 1 << ((ai * 4 + bj) * 4 + ci)
    return result


def schoolbook_scheme() -> MaskScheme:
    terms: list[MaskTerm] = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ai = i * 2 + k
                bj = k * 2 + j
                ci = i * 2 + j
                terms.append((1 << ai, 1 << bj, 1 << ci))
    result = canonical_scheme(terms)
    assert tensor_of_scheme(result) == target_matrix_multiplication_tensor()
    return result


def _permute(term: MaskTerm, permutation: tuple[int, int, int]) -> MaskTerm:
    return term[permutation[0]], term[permutation[1]], term[permutation[2]]


def _unpermute(term: MaskTerm, permutation: tuple[int, int, int]) -> MaskTerm:
    inverse = [0, 0, 0]
    for index, value in enumerate(permutation):
        inverse[value] = index
    return term[inverse[0]], term[inverse[1]], term[inverse[2]]


def exact_reductions(scheme: MaskScheme) -> tuple[MaskScheme, ...]:
    """All one-step rank reductions, each checked against the input tensor."""
    generated: list[MaskScheme] = []
    terms = list(scheme)
    original_tensor = tensor_of_scheme(scheme)
    for x in range(len(terms)):
        for y in range(x + 1, len(terms)):
            for permutation in PERMUTATIONS:
                a, b, c = _permute(terms[x], permutation)
                a2, b2, c2 = _permute(terms[y], permutation)
                if a != a2 or b != b2:
                    continue
                merged_c = c ^ c2
                remaining = [terms[i] for i in range(len(terms)) if i not in (x, y)]
                if merged_c:
                    remaining.append(_unpermute((a, b, merged_c), permutation))
                candidate = canonical_scheme(remaining)
                if len(candidate) >= len(scheme):
                    continue
                if tensor_of_scheme(candidate) != original_tensor:
                    raise AssertionError("invalid GF(2) rank reduction generated")
                generated.append(candidate)
    return tuple(dict.fromkeys(generated))


def exact_flips(scheme: MaskScheme) -> tuple[MaskScheme, ...]:
    """All exact, rank-preserving local flips adjacent to ``scheme``."""
    generated: list[MaskScheme] = []
    terms = list(scheme)
    original_tensor = tensor_of_scheme(scheme)
    for x in range(len(terms)):
        for y in range(x + 1, len(terms)):
            for permutation in PERMUTATIONS:
                a, b, c = _permute(terms[x], permutation)
                a2, b2, c2 = _permute(terms[y], permutation)
                if a != a2 or b == b2 or c == c2:
                    continue
                new_c = c ^ c2
                new_b = b2 ^ b
                if not new_c or not new_b:
                    continue
                first = _unpermute((a, b, new_c), permutation)
                second = _unpermute((a, new_b, c2), permutation)
                remaining = [terms[i] for i in range(len(terms)) if i not in (x, y)]
                remaining.extend((first, second))
                candidate = canonical_scheme(remaining)
                if len(candidate) != len(scheme):
                    continue
                if tensor_of_scheme(candidate) != original_tensor:
                    raise AssertionError("invalid GF(2) flip generated")
                generated.append(candidate)
    return tuple(dict.fromkeys(generated))


def naive_xor_addition_cost(scheme: MaskScheme) -> int:
    """Naive linear-form cost; used only as a deterministic tie-breaker."""
    input_cost = sum(
        max(0, a.bit_count() - 1) + max(0, b.bit_count() - 1)
        for a, b, _ in scheme
    )
    output_uses = [0] * 4
    for _, _, c in scheme:
        for output in range(4):
            if (c >> output) & 1:
                output_uses[output] += 1
    output_cost = sum(max(0, uses - 1) for uses in output_uses)
    return input_cost + output_cost


def search_rank7(seed: int = 0, max_flips: int = 200_000) -> FlipSearchResult:
    """Perform an exact random walk until a rank-7 decomposition is reached."""
    start_time = perf_counter()
    rng = random.Random(seed)
    scheme = schoolbook_scheme()
    target = target_matrix_multiplication_tensor()
    visited = {scheme}

    for step in range(max_flips + 1):
        reductions = exact_reductions(scheme)
        if reductions:
            candidate = min(reductions, key=naive_xor_addition_cost)
            if len(candidate) != 7 or tensor_of_scheme(candidate) != target:
                raise AssertionError("search returned an invalid rank-7 support")
            return FlipSearchResult(
                seed=seed,
                scheme=candidate,
                flips=step,
                visited_states=len(visited),
                elapsed_seconds=perf_counter() - start_time,
                naive_xor_additions=naive_xor_addition_cost(candidate),
            )

        if step == max_flips:
            break
        neighbors = exact_flips(scheme)
        if not neighbors:
            raise RuntimeError("flip graph walk reached a dead end")
        unseen = [candidate for candidate in neighbors if candidate not in visited]
        scheme = rng.choice(unseen if unseen else neighbors)
        visited.add(scheme)

    raise TimeoutError(f"no rank-7 scheme found within {max_flips} flips for seed {seed}")


def search_best(seeds: Iterable[int], max_flips: int = 200_000) -> tuple[FlipSearchResult, list[FlipSearchResult]]:
    results = [search_rank7(seed, max_flips) for seed in seeds]
    best = min(results, key=lambda item: (item.naive_xor_additions, item.flips, item.seed))
    return best, results
