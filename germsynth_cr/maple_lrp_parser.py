from __future__ import annotations

import bz2
import re
from fractions import Fraction
from pathlib import Path

from .exact_tensor import Scheme, SparseRow


MATRIX = re.compile(r"Matrix\((\d+),\s*(\d+),\s*\[\[")


def _read_matrix(source: str, offset: int):
    """Independent LRP tokenizer; intentionally does not import the Tensor parser."""
    header = MATRIX.match(source, offset)
    if not header:
        raise ValueError("LRP Matrix header missing")
    height, width = map(int, header.groups())
    cursor = header.end()
    rows = []
    for expected_row in range(height):
        row: SparseRow = {}
        for column in range(width):
            end = cursor
            while source[end] not in ",]":
                end += 1
            value = Fraction(source[cursor:end].strip())
            if value:
                row[column] = value
            delimiter = source[end]
            cursor = end + 1
            if column + 1 < width and delimiter != ",":
                raise ValueError("LRP row ended early")
        rows.append(row)
        if expected_row + 1 < height:
            if source[cursor:cursor + 2] != ",[":
                raise ValueError("LRP row separator missing")
            cursor += 2
    return height, width, rows


def parse_lrp(path: str | Path, shape: tuple[int, int, int]) -> Scheme:
    path = Path(path)
    text = bz2.open(path, "rt").read() if path.suffix == ".bz2" else path.read_text()
    starts = [match.start() for match in MATRIX.finditer(text)]
    if len(starts) != 3:
        raise ValueError(f"LRP expected three matrices, found {len(starts)}")
    matrices = [_read_matrix(text, start) for start in starts]
    m, n, p = shape
    rank = matrices[0][0]
    if tuple((x[0], x[1]) for x in matrices) != ((rank, m * n), (rank, n * p), (m * p, rank)):
        raise ValueError("LRP factor dimensions mismatch")
    U = [{(index % m) * n + index // m: value for index, value in row.items()} for row in matrices[0][2]]
    V = [{(index % n) * p + index // n: value for index, value in row.items()} for row in matrices[1][2]]
    W: list[SparseRow] = [{} for _ in range(rank)]
    for coordinate, row in enumerate(matrices[2][2]):
        for product, value in row.items():
            W[product][coordinate] = value
    scheme = Scheme("Q", shape, rank, U, V, W, str(path),
                    {"parser": "germsynth_cr.maple_lrp_parser", "representation": "LRP"})
    scheme.validate_dimensions()
    return scheme
