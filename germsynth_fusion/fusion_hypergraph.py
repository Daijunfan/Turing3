from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass


@dataclass
class HypergraphComponent:
    slot_ids: list[int]
    product_ids: list[int]


def connected_components(product_slot_incidence: dict[int, set[int]]) -> list[HypergraphComponent]:
    slot_products: dict[int, set[int]] = defaultdict(set)
    for product, slots in product_slot_incidence.items():
        for slot in slots:
            slot_products[slot].add(product)
    seen_slots: set[int] = set()
    seen_products: set[int] = set()
    result: list[HypergraphComponent] = []
    for seed in sorted(slot_products):
        if seed in seen_slots:
            continue
        slots: set[int] = set()
        products: set[int] = set()
        queue = deque([("slot", seed)])
        while queue:
            kind, node = queue.popleft()
            if kind == "slot":
                if node in seen_slots:
                    continue
                seen_slots.add(node)
                slots.add(node)
                queue.extend(("product", product) for product in slot_products[node])
            else:
                if node in seen_products:
                    continue
                seen_products.add(node)
                products.add(node)
                queue.extend(("slot", slot) for slot in product_slot_incidence[node])
        result.append(HypergraphComponent(sorted(slots), sorted(products)))
    return result


def canonical_component(component: HypergraphComponent) -> tuple[int, tuple[int, ...]]:
    """A product-id-independent coarse key used for symmetric census de-duplication."""
    incidence_sizes = tuple(sorted(len(component.slot_ids) for _ in component.product_ids))
    return len(component.slot_ids), incidence_sizes
