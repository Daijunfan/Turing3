from __future__ import annotations

import bz2
import re
from fractions import Fraction
from pathlib import Path

from .exact_tensor import Scheme, SparseRow


HEADER = re.compile(r"Matrix\((\d+),\s*(\d+),\s*\[\[")


def _matrix(text: str, start: int):
    match = HEADER.match(text, start)
    if not match:
        raise ValueError("Tensor parser expected Matrix")
    height, width = map(int, match.groups())
    pos, token_start, column = match.end(), match.end(), 0
    rows: list[SparseRow] = []
    row: SparseRow = {}
    while pos < len(text):
        char = text[pos]
        if char in ",]":
            token = text[token_start:pos].strip()
            if token:
                value = Fraction(token)
                if value:
                    row[column] = value
                column += 1
            if char == "]":
                if column != width:
                    raise ValueError(f"Tensor matrix row width {column} != {width}")
                rows.append(row)
                row, column = {}, 0
                if text.startswith("])", pos + 1):
                    if len(rows) != height:
                        raise ValueError("Tensor matrix height mismatch")
                    return height, width, rows
                if text.startswith(",[", pos + 1):
                    pos += 2
                    token_start = pos + 1
            else:
                token_start = pos + 1
        pos += 1
    raise ValueError("unterminated Tensor matrix")


def _flatten(rows: list[SparseRow], width: int) -> SparseRow:
    return {i * width + j: value for i, row in enumerate(rows) for j, value in row.items()}


def parse_tensor(path: str | Path, shape: tuple[int, int, int] | None = None) -> Scheme:
    path = Path(path)
    text = bz2.open(path, "rt").read() if path.suffix == ".bz2" else path.read_text()
    if shape is None:
        dims = re.findall(r"^[ABC]:=Matrix\((\d+),\s*(\d+),", text, re.MULTILINE)
        if len(dims) < 3:
            raise ValueError("Tensor source matrix declarations missing")
        shape = (int(dims[0][0]), int(dims[0][1]), int(dims[1][1]))
    marker = text.find("Tensor:=TriadSet([")
    if marker < 0:
        raise ValueError("Tensor TriadSet missing")
    body = text[marker:]
    starts = [match.start() for match in HEADER.finditer(body)]
    if len(starts) % 3:
        raise ValueError("incomplete Tensor triad")
    m, n, p = shape
    U: list[SparseRow] = []
    V: list[SparseRow] = []
    W: list[SparseRow] = []
    for product in range(len(starts) // 3):
        a, b, c = (_matrix(body, starts[3 * product + offset]) for offset in range(3))
        if (a[0], a[1], b[0], b[1], c[0], c[1]) != (m, n, n, p, p, m):
            raise ValueError(f"Tensor triad {product} dimensions mismatch")
        U.append(_flatten(a[2], n))
        V.append(_flatten(b[2], p))
        W.append({i * p + j: value for j, row in enumerate(c[2]) for i, value in row.items()})
    scheme = Scheme("Q", shape, len(U), U, V, W, str(path),
                    {"parser": "germsynth_cr.maple_tensor_parser", "representation": "TriadSet"})
    scheme.validate_dimensions()
    return scheme
