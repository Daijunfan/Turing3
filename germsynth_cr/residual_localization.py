from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction

from .residual import Residual


def sparse_rank(rows: list[dict[object, Fraction]]) -> int:
    pivots: dict[object, dict[object, Fraction]] = {}
    for source in rows:
        row = dict(source)
        while row:
            pivot = min(row)
            coefficient = row[pivot]
            if pivot not in pivots:
                pivots[pivot] = {column: value / coefficient for column, value in row.items()}
                break
            basis = pivots[pivot]
            for column, value in basis.items():
                updated = row.get(column, Fraction()) - coefficient * value
                if updated:
                    row[column] = updated
                else:
                    row.pop(column, None)
    return len(pivots)


def flattening_ranks(residual: Residual) -> tuple[int, int, int]:
    m, n, p = residual.shape
    sizes = (m * n, n * p, m * p)
    ranks = []
    for mode, size in enumerate(sizes):
        rows: list[dict[object, Fraction]] = [{} for _ in range(size)]
        for coordinate, value in residual.values.items():
            row = coordinate[mode]
            column = tuple(coordinate[index] for index in range(3) if index != mode)
            rows[row][column] = value
        ranks.append(sparse_rank(rows))
    return tuple(ranks)


def support_components(residual: Residual) -> list[dict]:
    parent = {}
    def find(node):
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node
    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a
    for a, b, c in residual.values:
        nodes = ((0, a), (1, b), (2, c))
        union(nodes[0], nodes[1]); union(nodes[0], nodes[2])
    grouped = defaultdict(lambda: {"coordinates": 0, "axis_support": [set(), set(), set()]})
    for coordinate in residual.values:
        root = find((0, coordinate[0]))
        grouped[root]["coordinates"] += 1
        for axis, index in enumerate(coordinate):
            grouped[root]["axis_support"][axis].add(index)
    return sorted(({"nonzero_count": value["coordinates"],
                    "axis_support": [sorted(support) for support in value["axis_support"]]}
                   for value in grouped.values()), key=lambda item: (-item["nonzero_count"], item["axis_support"]))


def localize_residual(residual: Residual) -> dict:
    ranks = flattening_ranks(residual)
    axis_support = [sorted({coordinate[axis] for coordinate in residual.values}) for axis in range(3)]
    return {
        "nonzero_count": len(residual.values),
        "coefficient_distribution": dict(sorted(Counter(str(value) for value in residual.values.values()).items())),
        "flattening_ranks": list(ranks),
        "flattening_lower_bound": max(ranks, default=0),
        "axis_support": axis_support,
        "axis_support_sizes": list(map(len, axis_support)),
        "connected_support_components": support_components(residual),
        "residual_sha256": residual.sha256(),
    }
