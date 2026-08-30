from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .fusion_hypergraph import HypergraphComponent, connected_components
from .kron_structure import column_space_basis, factor_kron_ranks, in_span, reshape_factor
from .scheme_io import Scheme
from .slot_inference import infer_slots


@dataclass
class Extraction:
    ordinary_slots: list[int]
    ordinary_product_ids: list[int]
    exceptional_slots: list[int]
    components: list[HypergraphComponent]
    product_slot_incidence: dict[int, set[int]]
    kron_ranks: list[tuple[int, int, int]]


def _direction(vector: list[Fraction]) -> tuple[Fraction, ...]:
    pivot = next((value for value in vector if value), None)
    return tuple(value / pivot for value in vector) if pivot is not None else ()


def _outer_intersection(parent: Scheme, outer: Scheme, inner: Scheme, product: int) -> set[int]:
    m, n, p = parent.shape
    im, inn, ip = inner.shape
    result: list[set[int]] = []
    for factor, rows, dims, width in zip(
        (parent.U[product], parent.V[product], parent.W[product]),
        (outer.U, outer.V, outer.W),
        ((m, n, im, inn), (n, p, inn, ip), (m, p, im, ip)),
        (outer.shape[0] * outer.shape[1], outer.shape[1] * outer.shape[2], outer.shape[0] * outer.shape[2]),
    ):
        basis = column_space_basis(reshape_factor(factor, *dims))
        result.append({slot for slot, row in enumerate(rows) if in_span(row, basis, width)})
    return set.intersection(*result)


def _factor_candidates(
    parent: Scheme, outer: Scheme, inner: Scheme, product: int, slots: set[int]
) -> tuple[set[int], set[int], set[int]]:
    m, n, p = parent.shape
    im, inn, ip = inner.shape
    result = []
    for factor, rows, dims, width in zip(
        (parent.U[product], parent.V[product], parent.W[product]),
        (outer.U, outer.V, outer.W),
        ((m, n, im, inn), (n, p, inn, ip), (m, p, im, ip)),
        (outer.shape[0] * outer.shape[1], outer.shape[1] * outer.shape[2], outer.shape[0] * outer.shape[2]),
    ):
        basis = column_space_basis(reshape_factor(factor, *dims))
        result.append({slot for slot in slots if in_span(rows[slot], basis, width)})
    return tuple(result)  # type: ignore[return-value]


def extract_components(parent: Scheme, outer: Scheme, inner: Scheme) -> Extraction:
    """Recover ordinary slots first, then close all remaining products over exceptional slots."""
    inferred = infer_slots(parent, outer, inner)
    ranks = inferred.kron_ranks
    owners: dict[int, int] = {}
    slot_counts: Counter[int] = Counter()
    for product, triple in enumerate(ranks):
        if triple != (1, 1, 1):
            continue
        candidates = set(inferred.product_slot_incidence[product])
        if len(candidates) == 1:
            owner = next(iter(candidates))
            owners[product] = owner
            slot_counts[owner] += 1
    ordinary_slots = sorted(slot for slot in range(outer.rank) if slot_counts[slot] == inner.rank)
    ordinary_set = set(ordinary_slots)
    ordinary_products = sorted(product for product, slot in owners.items() if slot in ordinary_set)
    exceptional_slots = set(range(outer.rank)) - ordinary_set
    exceptional_products = set(range(parent.rank)) - set(ordinary_products)
    incidence: dict[int, set[int]] = {}
    for product in exceptional_products:
        per_factor = _factor_candidates(parent, outer, inner, product, exceptional_slots)
        slots = set.union(*per_factor)
        if slots:
            incidence[product] = slots
    components = connected_components(incidence)
    assigned = {product for component in components for product in component.product_ids}
    # Products with no recognizable outer direction remain a separate whole-scheme-cancellation component.
    unassigned = sorted(exceptional_products - assigned)
    if unassigned:
        components.append(HypergraphComponent([], unassigned))
    return Extraction(ordinary_slots, ordinary_products, sorted(exceptional_slots), components, incidence, ranks)


def _add_product_tensor(
    accum: dict[int, dict[tuple[int, int], Fraction]], U, V, W, scale: Fraction = Fraction(1)
) -> None:
    for c, wc in W.items():
        out = accum[c]
        for a, uc in U.items():
            coefficient = scale * uc * wc
            for b, vc in V.items():
                key = (a, b)
                value = out.get(key, Fraction(0)) + coefficient * vc
                if value:
                    out[key] = value
                else:
                    out.pop(key, None)


def _embedded_slot_target(
    accum: dict[int, dict[tuple[int, int], Fraction]], outer: Scheme, inner: Scheme, slot: int
) -> None:
    om, on, op = outer.shape
    im, inn, ip = inner.shape
    pm, pn, pp = om * im, on * inn, op * ip
    for ao, uo in outer.U[slot].items():
        oi, ok = divmod(ao, on)
        for bo, vo in outer.V[slot].items():
            ok2, oj = divmod(bo, op)
            if ok != ok2:
                # Outer factors are arbitrary linear forms: indices need not themselves form a target triple.
                pass
            for co, wo in outer.W[slot].items():
                oi2, oj2 = divmod(co, op)
                outer_coefficient = uo * vo * wo
                for ii in range(im):
                    for ik in range(inn):
                        a = (oi * im + ii) * pn + (ok * inn + ik)
                        for ij in range(ip):
                            b = (ok2 * inn + ik) * pp + (oj * ip + ij)
                            c = (oi2 * im + ii) * pp + (oj2 * ip + ij)
                            out = accum[c]
                            key = (a, b)
                            value = out.get(key, Fraction(0)) - outer_coefficient
                            if value:
                                out[key] = value
                            else:
                                out.pop(key, None)


def component_residual(
    parent: Scheme, outer: Scheme, inner: Scheme, component: HypergraphComponent, sample_limit: int = 64
) -> dict[str, Any]:
    accum = component_residual_map(parent, outer, inner, component)
    residual_count = sum(1 for out in accum.values() for value in out.values() if value)
    sample = []
    for c in sorted(accum):
        for (a, b), value in sorted(accum[c].items()):
            if value and len(sample) < sample_limit:
                sample.append({"a": a, "b": b, "c": c, "value": str(value)})
    isolated_rank = len(component.slot_ids) * inner.rank
    return {
        "slot_ids": component.slot_ids,
        "product_ids": component.product_ids,
        "slot_count": len(component.slot_ids),
        "product_count": len(component.product_ids),
        "isolated_rank": isolated_rank,
        "fusion_gain": isolated_rank - len(component.product_ids),
        "residual_count": residual_count,
        "residual_sample": sample,
        "standalone": residual_count == 0 and bool(component.slot_ids),
        "classification": "standalone local kernel" if residual_count == 0 and component.slot_ids
                          else "whole-scheme cancellation" if not component.slot_ids
                          else "invalid or parent-dependent kernel",
    }


def component_residual_map(
    parent: Scheme, outer: Scheme, inner: Scheme, component: HypergraphComponent
) -> dict[int, dict[tuple[int, int], Fraction]]:
    accum: dict[int, dict[tuple[int, int], Fraction]] = defaultdict(dict)
    for product in component.product_ids:
        _add_product_tensor(accum, parent.U[product], parent.V[product], parent.W[product])
    for slot in component.slot_ids:
        _embedded_slot_target(accum, outer, inner, slot)
    return accum


def dependency_closures(
    parent: Scheme, outer: Scheme, inner: Scheme, components: list[HypergraphComponent]
) -> list[dict[str, Any]]:
    """Close components whose residual supports touch; report the exact merged residual."""
    maps = [component_residual_map(parent, outer, inner, component) for component in components]
    supports = [{(a, b, c) for c, out in residual.items() for (a, b), value in out.items() if value}
                for residual in maps]
    adjacency = [set() for _ in components]
    for i in range(len(components)):
        for j in range(i + 1, len(components)):
            if supports[i] & supports[j]:
                adjacency[i].add(j)
                adjacency[j].add(i)
    seen = set()
    result = []
    for seed in range(len(components)):
        if seed in seen:
            continue
        stack = [seed]
        group = []
        while stack:
            index = stack.pop()
            if index in seen:
                continue
            seen.add(index)
            group.append(index)
            stack.extend(adjacency[index] - seen)
        merged = HypergraphComponent(
            sorted({slot for index in group for slot in components[index].slot_ids}),
            sorted({product for index in group for product in components[index].product_ids}),
        )
        report = component_residual(parent, outer, inner, merged)
        report["source_component_indices"] = sorted(group)
        report["minimal_dependency_closure"] = True
        result.append(report)
    return result
