"""Exact selection of a small regenerative basis of algorithm phenotypes.

A phenotype requires a set of abstract local resources (here: linear-form
supports).  A pool survives a failure set if at least one phenotype avoids all
failed resources.  The routines below solve the resulting finite cover problem
exactly for small candidate pools.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Sequence

from .polymorphic import PhenotypePool, Resource, phenotype_resources


@dataclass(frozen=True, slots=True)
class CoverResult:
    phenotype_indices: tuple[int, ...]
    failure_order: int
    covered_failure_sets: int
    total_failure_sets: int

    @property
    def coverage_fraction(self) -> float:
        return self.covered_failure_sets / self.total_failure_sets if self.total_failure_sets else 1.0


class RegenerativeCoverOptimizer:
    def __init__(self, pool: PhenotypePool, universe: Sequence[Resource] | None = None) -> None:
        self.pool = pool
        self.universe = tuple(universe) if universe is not None else pool.resource_universe
        self.requirements = tuple(phenotype_resources(germ) for germ in pool.germs)
        self._failure_masks: dict[int, tuple[int, ...]] = {}

    def failure_masks(self, failure_order: int) -> tuple[int, ...]:
        if failure_order < 0:
            raise ValueError("failure_order must be nonnegative")
        cached = self._failure_masks.get(failure_order)
        if cached is not None:
            return cached
        masks: list[int] = []
        for failed in combinations(self.universe, failure_order):
            feasible_mask = 0
            for index, requirements in enumerate(self.requirements):
                if requirements.isdisjoint(failed):
                    feasible_mask |= 1 << index
            masks.append(feasible_mask)
        result = tuple(masks)
        self._failure_masks[failure_order] = result
        return result

    def evaluate(self, phenotype_indices: Iterable[int], failure_order: int) -> CoverResult:
        indices = tuple(sorted(set(phenotype_indices)))
        if any(index < 0 or index >= len(self.pool.germs) for index in indices):
            raise IndexError("phenotype index out of range")
        selected_mask = sum(1 << index for index in indices)
        failure_masks = self.failure_masks(failure_order)
        covered = sum(bool(selected_mask & feasible_mask) for feasible_mask in failure_masks)
        return CoverResult(indices, failure_order, covered, len(failure_masks))

    def best_exact_size(self, size: int, failure_order: int) -> CoverResult:
        if size <= 0 or size > len(self.pool.germs):
            raise ValueError("invalid subset size")
        best: CoverResult | None = None
        for indices in combinations(range(len(self.pool.germs)), size):
            candidate = self.evaluate(indices, failure_order)
            if best is None or (
                candidate.covered_failure_sets,
                tuple(-index for index in candidate.phenotype_indices),
            ) > (
                best.covered_failure_sets,
                tuple(-index for index in best.phenotype_indices),
            ):
                best = candidate
        assert best is not None
        return best

    def minimum_full_cover(self, failure_order: int) -> CoverResult | None:
        total = len(self.failure_masks(failure_order))
        for size in range(1, len(self.pool.germs) + 1):
            best = self.best_exact_size(size, failure_order)
            if best.covered_failure_sets == total:
                return best
        return None

    def minimum_single_fault_cover_with_tiebreak(self) -> CoverResult:
        """Minimize pool size, then maximize pair and triple fault coverage."""
        single_masks = self.failure_masks(1)
        total_single = len(single_masks)
        for size in range(1, len(self.pool.germs) + 1):
            full_single: list[tuple[int, ...]] = []
            for indices in combinations(range(len(self.pool.germs)), size):
                result = self.evaluate(indices, 1)
                if result.covered_failure_sets == total_single:
                    full_single.append(indices)
            if not full_single:
                continue
            best_indices = max(
                full_single,
                key=lambda indices: (
                    self.evaluate(indices, 2).covered_failure_sets,
                    self.evaluate(indices, 3).covered_failure_sets,
                    tuple(-index for index in indices),
                ),
            )
            return self.evaluate(best_indices, 1)
        raise RuntimeError("no single-fault cover exists")
