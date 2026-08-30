from __future__ import annotations

from dataclasses import dataclass

from .kron_structure import column_space_basis, factor_kron_ranks, in_span, reshape_factor
from .scheme_io import Scheme


@dataclass
class SlotInference:
    kron_ranks: list[tuple[int, int, int]]
    product_slot_incidence: list[list[int]]
    mixed_product_ids: list[int]
    plain_product_ids: list[int]
    unresolved_product_ids: list[int]


def _direction_dense(vector) -> tuple:
    pivot = next((value for value in vector if value), None)
    return tuple(value / pivot for value in vector) if pivot is not None else ()


def _direction_sparse(row, width: int) -> tuple:
    return _direction_dense([row.get(i, 0) for i in range(width)])


def _candidate_slots(basis, direction_map, outer_rows, width):
    if len(basis) == 1:
        return direction_map.get(_direction_dense(basis[0]), set())
    return {slot for slot, row in outer_rows if in_span(row, basis, width)}


def infer_slots(parent: Scheme, outer: Scheme, inner: Scheme) -> SlotInference:
    expected = tuple(a * b for a, b in zip(outer.shape, inner.shape))
    if parent.shape != expected:
        raise ValueError(f"parent shape {parent.shape} != outer*inner {expected}")
    m, n, p = parent.shape
    im, inn, ip = inner.shape
    ranks = factor_kron_ranks(parent, inner.shape)
    outer_widths = (outer.shape[0] * outer.shape[1], outer.shape[1] * outer.shape[2],
                    outer.shape[0] * outer.shape[2])
    direction_maps = []
    for rows, width in zip((outer.U, outer.V, outer.W), outer_widths):
        mapping = {}
        for slot, row in enumerate(rows):
            mapping.setdefault(_direction_sparse(row, width), set()).add(slot)
        direction_maps.append(mapping)
    incidence: list[list[int]] = []
    mixed: list[int] = []
    plain: list[int] = []
    unresolved: list[int] = []
    for product, (U, V, W) in enumerate(zip(parent.U, parent.V, parent.W)):
        ub = column_space_basis(reshape_factor(U, m, n, im, inn))
        vb = column_space_basis(reshape_factor(V, n, p, inn, ip))
        wb = column_space_basis(reshape_factor(W, m, p, im, ip))
        candidates = [
            _candidate_slots(basis, direction_map, list(enumerate(rows)), width)
            for basis, direction_map, rows, width in zip(
                (ub, vb, wb), direction_maps, (outer.U, outer.V, outer.W), outer_widths)
        ]
        slots = sorted(set.intersection(*candidates))
        incidence.append(slots)
        if len(slots) >= 2:
            mixed.append(product)
        elif len(slots) == 1 and ranks[product] == (1, 1, 1):
            plain.append(product)
        else:
            unresolved.append(product)
    return SlotInference(ranks, incidence, mixed, plain, unresolved)
