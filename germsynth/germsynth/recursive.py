"""Arbitrary-scale algorithms generated from exact bilinear germs."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .germ import BilinearGerm

Matrix = list[list[int]]
Polynomial = list[int]


@dataclass(slots=True)
class OperationCounter:
    scalar_multiplications: int = 0
    ring_additions: int = 0
    negations: int = 0
    recursive_calls: int = 0


def _is_power_of_two(value: int) -> bool:
    return value > 0 and (value & (value - 1)) == 0


def _validate_square_matrix(matrix: Sequence[Sequence[int]]) -> int:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    if not _is_power_of_two(n):
        raise ValueError("matrix dimension must be a power of two")
    return n


def _zero_matrix(n: int) -> Matrix:
    return [[0 for _ in range(n)] for _ in range(n)]


def _split_matrix(matrix: Sequence[Sequence[int]]) -> tuple[Matrix, Matrix, Matrix, Matrix]:
    n = len(matrix)
    half = n // 2
    return (
        [list(row[:half]) for row in matrix[:half]],
        [list(row[half:]) for row in matrix[:half]],
        [list(row[:half]) for row in matrix[half:]],
        [list(row[half:]) for row in matrix[half:]],
    )


def _join_matrix(blocks: Sequence[Matrix]) -> Matrix:
    if len(blocks) != 4:
        raise ValueError("2x2 block matrix requires four blocks")
    top_left, top_right, bottom_left, bottom_right = blocks
    n = len(top_left)
    return [
        top_left[row] + top_right[row] for row in range(n)
    ] + [
        bottom_left[row] + bottom_right[row] for row in range(n)
    ]


def _linear_matrix_form(
    coefficients: Sequence[int], blocks: Sequence[Matrix], counter: OperationCounter
) -> Matrix:
    selected = [(coefficient, block) for coefficient, block in zip(coefficients, blocks, strict=True) if coefficient]
    if not selected:
        return _zero_matrix(len(blocks[0]))
    n = len(selected[0][1])
    first_coefficient, first_block = selected[0]
    result = [
        [first_coefficient * first_block[row][column] for column in range(n)]
        for row in range(n)
    ]
    if first_coefficient < 0:
        counter.negations += n * n
    for coefficient, block in selected[1:]:
        if coefficient < 0:
            counter.negations += n * n
        for row in range(n):
            for column in range(n):
                result[row][column] += coefficient * block[row][column]
                counter.ring_additions += 1
    return result


def _accumulate_matrix(
    destination: Matrix | None,
    source: Matrix,
    coefficient: int,
    counter: OperationCounter,
) -> Matrix:
    n = len(source)
    if destination is None:
        result = [
            [coefficient * source[row][column] for column in range(n)]
            for row in range(n)
        ]
        if coefficient < 0:
            counter.negations += n * n
        return result
    if coefficient < 0:
        counter.negations += n * n
    for row in range(n):
        for column in range(n):
            destination[row][column] += coefficient * source[row][column]
            counter.ring_additions += 1
    return destination


def germ_matrix_multiply(
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
    germ: BilinearGerm,
    counter: OperationCounter | None = None,
) -> Matrix:
    """Multiply square power-of-two matrices using a discovered 2x2 germ."""
    if germ.specification.kind != "matrix_multiplication" or germ.block_factor != 2:
        raise ValueError("germ is not a 2x2 matrix multiplication germ")
    n = _validate_square_matrix(left)
    if _validate_square_matrix(right) != n:
        raise ValueError("matrix dimensions differ")
    if counter is None:
        counter = OperationCounter()
    counter.recursive_calls += 1
    if n == 1:
        counter.scalar_multiplications += 1
        return [[left[0][0] * right[0][0]]]

    left_blocks = _split_matrix(left)
    right_blocks = _split_matrix(right)
    products: list[Matrix] = []
    for term in germ.terms:
        left_form = _linear_matrix_form(term.u, left_blocks, counter)
        right_form = _linear_matrix_form(term.v, right_blocks, counter)
        products.append(germ_matrix_multiply(left_form, right_form, germ, counter))

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


def naive_matrix_multiply(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> Matrix:
    n = len(left)
    if n == 0 or len(right) != n or any(len(row) != n for row in left) or any(len(row) != n for row in right):
        raise ValueError("naive matrix multiply expects equally-sized square matrices")
    return [
        [sum(left[i][k] * right[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def _linear_polynomial_form(
    coefficients: Sequence[int], blocks: Sequence[Polynomial], counter: OperationCounter
) -> Polynomial:
    selected = [(coefficient, block) for coefficient, block in zip(coefficients, blocks, strict=True) if coefficient]
    if not selected:
        return [0] * len(blocks[0])
    first_coefficient, first_block = selected[0]
    result = [first_coefficient * value for value in first_block]
    if first_coefficient < 0:
        counter.negations += len(result)
    for coefficient, block in selected[1:]:
        if coefficient < 0:
            counter.negations += len(block)
        for index, value in enumerate(block):
            result[index] += coefficient * value
            counter.ring_additions += 1
    return result


def germ_polynomial_multiply(
    left: Sequence[int],
    right: Sequence[int],
    germ: BilinearGerm,
    counter: OperationCounter | None = None,
) -> Polynomial:
    """Multiply equal-length power-of-two coefficient arrays using a rank-3 germ."""
    if germ.specification.kind != "polynomial_multiplication" or germ.block_factor != 2:
        raise ValueError("germ is not a two-way polynomial multiplication germ")
    n = len(left)
    if n != len(right) or not _is_power_of_two(n):
        raise ValueError("polynomial lengths must be equal powers of two")
    if counter is None:
        counter = OperationCounter()
    counter.recursive_calls += 1
    if n == 1:
        counter.scalar_multiplications += 1
        return [left[0] * right[0]]

    half = n // 2
    left_blocks = (list(left[:half]), list(left[half:]))
    right_blocks = (list(right[:half]), list(right[half:]))
    products: list[Polynomial] = []
    for term in germ.terms:
        left_form = _linear_polynomial_form(term.u, left_blocks, counter)
        right_form = _linear_polynomial_form(term.v, right_blocks, counter)
        products.append(germ_polynomial_multiply(left_form, right_form, germ, counter))

    result = [0] * (2 * n - 1)
    touched = [False] * len(result)
    for term, product in zip(germ.terms, products, strict=True):
        for output_block, coefficient in enumerate(term.w):
            if not coefficient:
                continue
            offset = output_block * half
            if coefficient < 0:
                counter.negations += len(product)
            for index, value in enumerate(product):
                destination = offset + index
                if touched[destination]:
                    counter.ring_additions += 1
                result[destination] += coefficient * value
                touched[destination] = True
    return result


def naive_polynomial_multiply(left: Sequence[int], right: Sequence[int]) -> Polynomial:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result
