from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, getcontext
from math import log


getcontext().prec = 70


@dataclass(frozen=True)
class Edge:
    multiplicity: int
    source_type: int
    target_type: int
    scale: Decimal
    shape_transformation: str = "uniform"
    additive_work: int = 0


def matrix_at(edges: list[Edge], type_count: int, gamma: Decimal) -> list[list[Decimal]]:
    matrix = [[Decimal(0) for _ in range(type_count)] for _ in range(type_count)]
    for edge in edges:
        matrix[edge.source_type][edge.target_type] += Decimal(edge.multiplicity) * (edge.scale ** (-gamma))
    return matrix


def spectral_radius(matrix: list[list[Decimal]], iterations: int = 500) -> Decimal:
    n = len(matrix)
    # A+I is aperiodic on every recurrent class and has rho(A+I)=rho(A)+1.
    shifted = [[matrix[i][j] + (Decimal(1) if i == j else Decimal(0))
                for j in range(n)] for i in range(n)]
    vector = [Decimal(1) for _ in range(n)]
    eigenvalue = Decimal(0)
    for _ in range(iterations):
        nxt = [sum(shifted[i][j] * vector[j] for j in range(n)) for i in range(n)]
        norm = max(abs(value) for value in nxt)
        if not norm:
            return Decimal(0)
        nxt = [value / norm for value in nxt]
        if abs(norm - eigenvalue) < Decimal("1e-55"):
            return norm - 1
        vector, eigenvalue = nxt, norm
    return eigenvalue - 1


def solve_critical_exponent(edges: list[Edge], type_count: int, tolerance=Decimal("1e-30")) -> dict:
    low, high = Decimal(0), Decimal(10)
    while spectral_radius(matrix_at(edges, type_count, high)) > 1:
        high *= 2
    while high - low > tolerance:
        middle = (low + high) / 2
        if spectral_radius(matrix_at(edges, type_count, middle)) > 1:
            low = middle
        else:
            high = middle
    gamma = (low + high) / 2
    return {"gamma": str(gamma), "lower_bound": str(low), "upper_bound": str(high),
            "rho_lower": str(spectral_radius(matrix_at(edges, type_count, low))),
            "rho_upper": str(spectral_radius(matrix_at(edges, type_count, high)))}


def single_type_check(rank: int, scale: int) -> dict:
    solved = solve_critical_exponent([Edge(rank, 0, 0, Decimal(scale))], 1)
    expected = log(rank, scale)
    actual = float(solved["gamma"])
    return {**solved, "rank": rank, "scale": scale, "expected_log": expected,
            "absolute_error": abs(actual - expected), "pass": abs(actual - expected) <= 1e-12}


def empirical_single_type(rank: int, scale: int, levels: int = 12) -> dict:
    sizes = [scale ** level for level in range(1, levels + 1)]
    operations = [rank ** level for level in range(1, levels + 1)]
    xs = [log(size) for size in sizes]
    ys = [log(work) for work in operations]
    xbar, ybar = sum(xs) / levels, sum(ys) / levels
    gamma = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / sum((x - xbar) ** 2 for x in xs)
    theory = log(rank, scale)
    return {"levels": levels, "gamma_empirical": gamma, "gamma_theory": theory,
            "absolute_error": abs(gamma - theory), "pass": abs(gamma - theory) <= 0.01}
