from __future__ import annotations

from fractions import Fraction

from .scheme_io import SparseRow


def reshape_factor(
    row: SparseRow, parent_rows: int, parent_cols: int, inner_rows: int, inner_cols: int
) -> list[list[Fraction]]:
    if parent_rows % inner_rows or parent_cols % inner_cols:
        raise ValueError("inner shape does not divide parent factor shape")
    outer_rows, outer_cols = parent_rows // inner_rows, parent_cols // inner_cols
    matrix = [[Fraction(0) for _ in range(inner_rows * inner_cols)]
              for _ in range(outer_rows * outer_cols)]
    for flat, value in row.items():
        i, j = divmod(flat, parent_cols)
        outer = (i // inner_rows) * outer_cols + j // inner_cols
        inner = (i % inner_rows) * inner_cols + j % inner_cols
        matrix[outer][inner] = value
    return matrix


def exact_rank(matrix: list[list[Fraction]]) -> int:
    if not matrix:
        return 0
    work = [row[:] for row in matrix if any(row)]
    rank = 0
    column_count = len(matrix[0])
    for column in range(column_count):
        pivot = next((r for r in range(rank, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pv = work[rank][column]
        work[rank] = [value / pv for value in work[rank]]
        for r in range(len(work)):
            if r != rank and work[r][column]:
                scale = work[r][column]
                work[r] = [a - scale * b for a, b in zip(work[r], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def column_space_basis(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    if not matrix:
        return []
    columns = [[matrix[r][c] for r in range(len(matrix))] for c in range(len(matrix[0]))]
    basis: list[list[Fraction]] = []
    current_rank = 0
    for column in columns:
        candidate = basis + [column]
        rank = exact_rank([list(row) for row in zip(*candidate)])
        if rank > current_rank:
            basis.append(column)
            current_rank = rank
    return basis


def in_span(vector: SparseRow, basis: list[list[Fraction]], width: int) -> bool:
    dense = [vector.get(i, Fraction(0)) for i in range(width)]
    if not basis:
        return not any(dense)
    before = exact_rank([list(row) for row in zip(*basis)])
    after = exact_rank([list(row) for row in zip(*(basis + [dense]))])
    return before == after


def factor_kron_ranks(parent, inner_shape: tuple[int, int, int]) -> list[tuple[int, int, int]]:
    m, n, p = parent.shape
    im, inn, ip = inner_shape
    ranks = []
    for U, V, W in zip(parent.U, parent.V, parent.W):
        ranks.append((
            exact_rank(reshape_factor(U, m, n, im, inn)),
            exact_rank(reshape_factor(V, n, p, inn, ip)),
            exact_rank(reshape_factor(W, m, p, im, ip)),
        ))
    return ranks
