"""Exact enumeration that discovers a rank-3 polynomial multiplication germ."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from time import perf_counter

from .tensor import RankOneTerm, degree1_polynomial_multiplication_specification, outer


@dataclass(frozen=True, slots=True)
class KaratsubaSearchResult:
    terms: tuple[RankOneTerm, ...]
    distinct_rank_one_tensors: int
    examined_pairs: int
    elapsed_seconds: float


def _nonzero_ternary_forms(dimension: int) -> list[tuple[int, ...]]:
    return [
        tuple(values)
        for values in product((-1, 0, 1), repeat=dimension)
        if any(values)
    ]


def _canonical_sign(
    u: tuple[int, ...], v: tuple[int, ...], w: tuple[int, ...]
) -> RankOneTerm:
    u_list, v_list, w_list = list(u), list(v), list(w)
    if next(value for value in u_list if value) < 0:
        u_list = [-value for value in u_list]
        w_list = [-value for value in w_list]
    if next(value for value in v_list if value) < 0:
        v_list = [-value for value in v_list]
        w_list = [-value for value in w_list]
    return RankOneTerm(tuple(u_list), tuple(v_list), tuple(w_list))


def search_rank3() -> KaratsubaSearchResult:
    start_time = perf_counter()
    specification = degree1_polynomial_multiplication_specification()
    unique: dict[tuple[int, ...], RankOneTerm] = {}
    for u in _nonzero_ternary_forms(2):
        for v in _nonzero_ternary_forms(2):
            for w in _nonzero_ternary_forms(3):
                term = _canonical_sign(u, v, w)
                unique.setdefault(outer(term.u, term.v, term.w), term)

    candidates = list(unique.items())
    lookup = unique
    target = specification.tensor
    examined = 0
    for first_index, (first_tensor, first_term) in enumerate(candidates):
        for second_tensor, second_term in candidates[first_index:]:
            examined += 1
            needed = tuple(
                target[index] - first_tensor[index] - second_tensor[index]
                for index in range(len(target))
            )
            third_term = lookup.get(needed)
            if third_term is None:
                continue
            terms = (first_term, second_term, third_term)
            reconstructed = [0] * len(target)
            for term in terms:
                for index, value in enumerate(term.tensor()):
                    reconstructed[index] += value
            if tuple(reconstructed) != target:
                raise AssertionError("rank-3 search produced an invalid identity")
            return KaratsubaSearchResult(
                terms=terms,
                distinct_rank_one_tensors=len(candidates),
                examined_pairs=examined,
                elapsed_seconds=perf_counter() - start_time,
            )
    raise RuntimeError("no ternary rank-3 decomposition found")
