"""Locally reconfigurable exact algorithms built from equivalent germ phenotypes.

Every phenotype proves the same local tensor identity.  Therefore a recursive
node may choose any phenotype independently of every other node; correctness is
closed under arbitrary node-wise substitution.  This is the prototype's
computational analogue of functional degeneracy and local regeneration.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import random
from typing import Callable, Iterable, Sequence

from .germ import BilinearGerm
from .recursive import (
    Matrix,
    OperationCounter,
    _accumulate_matrix,
    _join_matrix,
    _linear_matrix_form,
    _split_matrix,
    _validate_square_matrix,
    _zero_matrix,
)

Resource = tuple[str, int]
Selector = Callable[[int, tuple[int, ...], int], int]


def _mix64(value: int) -> int:
    """Deterministic 64-bit mixer; independent of Python hash randomization."""
    mask = (1 << 64) - 1
    value = (value + 0x9E3779B97F4A7C15) & mask
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & mask
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & mask
    return value ^ (value >> 31)


def _stable_node_seed(seed: int, depth: int, path: tuple[int, ...], size: int) -> int:
    value = _mix64(seed & ((1 << 64) - 1))
    value = _mix64(value ^ depth)
    value = _mix64(value ^ size)
    for branch in path:
        value = _mix64(value ^ (branch + 1))
    return value


def _support_mask(coefficients: Sequence[int]) -> int:
    mask = 0
    for index, coefficient in enumerate(coefficients):
        if coefficient:
            mask |= 1 << index
    return mask


def phenotype_resources(germ: BilinearGerm) -> frozenset[Resource]:
    resources: set[Resource] = set()
    for term in germ.terms:
        resources.add(("left", _support_mask(term.u)))
        resources.add(("right", _support_mask(term.v)))
        resources.add(("output", _support_mask(term.w)))
    return frozenset(resources)


@dataclass(frozen=True, slots=True)
class PhenotypePool:
    germs: tuple[BilinearGerm, ...]

    def __post_init__(self) -> None:
        if not self.germs:
            raise ValueError("phenotype pool cannot be empty")
        first = self.germs[0]
        for germ in self.germs:
            if not germ.verify_local_identity():
                raise ValueError(f"unverified phenotype: {germ.name}")
            if germ.specification.tensor != first.specification.tensor:
                raise ValueError("phenotypes implement different tensors")
            if germ.block_factor != first.block_factor or germ.rank != first.rank:
                raise ValueError("phenotypes have different recursion shapes")

    @property
    def resource_universe(self) -> tuple[Resource, ...]:
        return tuple(sorted(set().union(*(phenotype_resources(germ) for germ in self.germs))))

    def feasible_indices(self, forbidden: Iterable[Resource]) -> tuple[int, ...]:
        blocked = frozenset(forbidden)
        return tuple(
            index
            for index, germ in enumerate(self.germs)
            if phenotype_resources(germ).isdisjoint(blocked)
        )

    def fault_coverage(self, failure_order: int = 1) -> dict[str, float | int]:
        universe = self.resource_universe
        total = 0
        fixed_survives = 0
        pool_survives = 0
        first_resources = phenotype_resources(self.germs[0])
        for failed in combinations(universe, failure_order):
            total += 1
            blocked = frozenset(failed)
            if first_resources.isdisjoint(blocked):
                fixed_survives += 1
            if self.feasible_indices(blocked):
                pool_survives += 1
        return {
            "failure_order": failure_order,
            "resource_universe": len(universe),
            "failure_sets": total,
            "fixed_survives": fixed_survives,
            "pool_survives": pool_survives,
            "fixed_survival_fraction": fixed_survives / total if total else 1.0,
            "pool_survival_fraction": pool_survives / total if total else 1.0,
        }


def random_selector(phenotype_count: int, seed: int) -> Selector:
    rng = random.Random(seed)

    def choose(_depth: int, _path: tuple[int, ...], _size: int) -> int:
        return rng.randrange(phenotype_count)

    return choose


def fault_aware_selector(
    pool: PhenotypePool,
    seed: int,
    failures_per_node: int = 1,
    *,
    strict: bool = False,
    failure_universe: Sequence[Resource] | None = None,
) -> Selector:
    """Choose a phenotype after deterministic node-local resource failures.

    Failures are sampled from an abstract linear-form resource universe.  With
    ``strict=True`` an infeasible sampled set is an error; otherwise progressively
    smaller failure sets are tried.  This models local reconfiguration, not
    physical gate timing.
    """
    if failures_per_node < 0:
        raise ValueError("failures_per_node must be nonnegative")
    universe = tuple(failure_universe) if failure_universe is not None else pool.resource_universe

    def choose(depth: int, path: tuple[int, ...], size: int) -> int:
        # Node identity creates a reproducible local stream independent of traversal order
        # and of Python's process-level hash seed.
        rng = random.Random(_stable_node_seed(seed, depth, path, size))
        requested_order = min(failures_per_node, len(universe))
        blocked = rng.sample(universe, requested_order) if requested_order else []
        feasible = pool.feasible_indices(blocked)
        if feasible:
            return rng.choice(feasible)
        if strict:
            raise RuntimeError(f"no phenotype survives local failures: {blocked}")
        for order in range(requested_order - 1, -1, -1):
            blocked = rng.sample(universe, order) if order else []
            feasible = pool.feasible_indices(blocked)
            if feasible:
                return rng.choice(feasible)
        raise RuntimeError("no feasible phenotype")

    return choose


def polymorphic_matrix_multiply(
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
    pool: PhenotypePool,
    selector: Selector,
    counter: OperationCounter | None = None,
    *,
    _depth: int = 0,
    _path: tuple[int, ...] = (),
) -> Matrix:
    n = _validate_square_matrix(left)
    if _validate_square_matrix(right) != n:
        raise ValueError("matrix dimensions differ")
    if counter is None:
        counter = OperationCounter()
    counter.recursive_calls += 1
    if n == 1:
        counter.scalar_multiplications += 1
        return [[left[0][0] * right[0][0]]]

    phenotype_index = selector(_depth, _path, n)
    if phenotype_index < 0 or phenotype_index >= len(pool.germs):
        raise IndexError("selector returned an invalid phenotype index")
    germ = pool.germs[phenotype_index]
    left_blocks = _split_matrix(left)
    right_blocks = _split_matrix(right)
    products: list[Matrix] = []
    for term_index, term in enumerate(germ.terms):
        left_form = _linear_matrix_form(term.u, left_blocks, counter)
        right_form = _linear_matrix_form(term.v, right_blocks, counter)
        products.append(
            polymorphic_matrix_multiply(
                left_form,
                right_form,
                pool,
                selector,
                counter,
                _depth=_depth + 1,
                _path=_path + (term_index,),
            )
        )

    output_blocks: list[Matrix | None] = [None, None, None, None]
    for term, product in zip(germ.terms, products, strict=True):
        for output_index, coefficient in enumerate(term.w):
            if coefficient:
                output_blocks[output_index] = _accumulate_matrix(
                    output_blocks[output_index], product, coefficient, counter
                )
    block_size = n // 2
    finalized = [block if block is not None else _zero_matrix(block_size) for block in output_blocks]
    return _join_matrix(finalized)
