from __future__ import annotations

from dataclasses import dataclass

from .contact_rank import ContactOrgan


@dataclass
class ContactHypergraph:
    slot_count: int
    edges: list[ContactOrgan]


class ContactPackingSolver:
    """Exact dynamic programming over disjoint hyperedges; shared slots are rejected."""
    def solve(self, graph: ContactHypergraph):
        by_first = {slot: [] for slot in range(graph.slot_count)}
        for edge in graph.edges:
            if not edge.exact:
                continue
            for slot in edge.slot_ids:
                by_first[slot].append(edge)
        memo = {}
        def visit(mask: int):
            if mask == (1 << graph.slot_count) - 1:
                return 0, []
            if mask in memo:
                return memo[mask]
            slot = next(index for index in range(graph.slot_count) if not (mask >> index) & 1)
            best = None
            for edge in by_first[slot]:
                edge_mask = sum(1 << index for index in edge.slot_ids)
                if edge_mask & mask:
                    continue
                tail = visit(mask | edge_mask)
                candidate = (edge.rank + tail[0], [edge] + tail[1])
                if best is None or candidate[0] < best[0]:
                    best = candidate
            if best is None:
                raise ValueError(f"no exact-cover edge for slot {slot}")
            memo[mask] = best
            return best
        cost, edges = visit(0)
        covered = [slot for edge in edges for slot in edge.slot_ids]
        if len(covered) != len(set(covered)) or set(covered) != set(range(graph.slot_count)):
            raise AssertionError("packing is not an exact cover")
        return {"cost": cost, "edges": edges}
