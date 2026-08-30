from __future__ import annotations

import bz2
import hashlib
import json
import re
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


Scalar = Fraction
SparseRow = dict[int, Scalar]


def _fraction(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        return Fraction(value)
    raise TypeError(f"unsupported exact coefficient: {value!r}")


def _sparse_row(values: Iterable[Any], modulus: int | None = None) -> SparseRow:
    row: SparseRow = {}
    for index, raw in enumerate(values):
        value = _fraction(raw)
        if modulus is not None:
            if value.denominator % modulus == 0:
                raise ValueError(f"denominator {value.denominator} is zero in GF({modulus})")
            value = Fraction((value.numerator * pow(value.denominator, -1, modulus)) % modulus)
        if value:
            row[index] = value
    return row


def _json_row(row: SparseRow) -> list[list[int | str]]:
    result: list[list[int | str]] = []
    for index, value in sorted(row.items()):
        encoded: int | str = value.numerator if value.denominator == 1 else str(value)
        result.append([index, encoded])
    return result


@dataclass
class Scheme:
    field: str
    shape: tuple[int, int, int]
    rank: int
    U: list[SparseRow]
    V: list[SparseRow]
    W: list[SparseRow]
    source: str
    provenance: dict[str, Any] = field(default_factory=dict)
    exact_coefficient_type: str = "rational"
    vectorization: dict[str, str] = field(default_factory=lambda: {
        "U": "row-major A[m,n]",
        "V": "row-major B[n,p]",
        "W": "row-major C[m,p]",
    })

    def validate_dimensions(self) -> None:
        m, n, p = self.shape
        if not (len(self.U) == len(self.V) == len(self.W) == self.rank):
            raise ValueError("declared rank does not match factor row counts")
        for name, rows, width in (("U", self.U, m * n), ("V", self.V, n * p), ("W", self.W, m * p)):
            for product, row in enumerate(rows):
                if any(index < 0 or index >= width for index in row):
                    raise ValueError(f"{name}[{product}] has coordinate outside width {width}")

    def to_json(self) -> dict[str, Any]:
        return {
            "format": "germsynth-fusion-scheme-v1",
            "field": self.field,
            "shape": list(self.shape),
            "rank": self.rank,
            "U": [_json_row(row) for row in self.U],
            "V": [_json_row(row) for row in self.V],
            "W": [_json_row(row) for row in self.W],
            "source": self.source,
            "provenance": self.provenance,
            "exact_coefficient_type": self.exact_coefficient_type,
            "vectorization": self.vectorization,
        }


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_sparse_json_rows(raw: list[Any], width: int, modulus: int | None) -> list[SparseRow]:
    rows: list[SparseRow] = []
    for item in raw:
        if isinstance(item, dict):
            dense = [0] * width
            for index, value in item.items():
                dense[int(index)] = value
            row = _sparse_row(dense, modulus)
        elif item and isinstance(item[0], list) and len(item[0]) == 2:
            dense = [0] * width
            for index, value in item:
                dense[int(index)] = value
            row = _sparse_row(dense, modulus)
        else:
            if len(item) != width:
                raise ValueError(f"dense factor row has width {len(item)}, expected {width}")
            row = _sparse_row(item, modulus)
        rows.append(row)
    return rows


def load_json_scheme(path: str | Path) -> Scheme:
    path = Path(path)
    raw = json.loads(path.read_text())
    shape = tuple(raw.get("shape", raw.get("n")))
    if len(shape) != 3:
        raise ValueError("scheme shape must have three dimensions")
    rank = int(raw.get("rank", raw.get("m")))
    field_name = raw.get("field") or ("GF(2)" if raw.get("z2") else "Q")
    modulus = 2 if field_name in ("GF(2)", "F2") else 3 if field_name in ("GF(3)", "F3") else None
    m, n, p = shape
    U = _parse_sparse_json_rows(raw.get("U", raw.get("u")), m * n, modulus)
    V = _parse_sparse_json_rows(raw.get("V", raw.get("v")), n * p, modulus)
    W = _parse_sparse_json_rows(raw.get("W", raw.get("w")), m * p, modulus)
    # Perminov stores W over C^T: transpose it into the unified row-major C convention.
    if "w" in raw and "W" not in raw:
        W = [{(index % m) * p + index // m: value for index, value in row.items()} for row in W]
    scheme = Scheme(
        field=str(field_name), shape=shape, rank=rank, U=U, V=V, W=W,
        source=str(raw.get("source", path)), provenance=raw.get("provenance", {"sha256": sha256_file(path)}),
        exact_coefficient_type=raw.get("exact_coefficient_type", "finite-field" if modulus else "rational"),
        vectorization=raw.get("vectorization", {
            "U": "row-major A[m,n]", "V": "row-major B[n,p]", "W": "row-major C[m,p]"
        }),
    )
    scheme.validate_dimensions()
    return scheme


_MATRIX_HEADER = re.compile(r"Matrix\((\d+),\s*(\d+),\s*\[\[")


def _parse_maple_matrix(text: str, start: int) -> tuple[int, int, list[SparseRow], int]:
    match = _MATRIX_HEADER.match(text, start)
    if not match:
        raise ValueError(f"expected Maple Matrix at offset {start}")
    height, width = map(int, match.groups())
    pos = match.end()
    rows: list[SparseRow] = []
    row: SparseRow = {}
    column = 0
    token_start = pos
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
                    raise ValueError(f"Maple row {len(rows)} has width {column}, expected {width}")
                rows.append(row)
                row = {}
                column = 0
                # Rows end with ], next row starts with ,[; matrix ends with ]]).
                if text.startswith("])", pos + 1):
                    end = pos + 3
                    if len(rows) != height:
                        raise ValueError(f"Maple matrix has {len(rows)} rows, expected {height}")
                    return height, width, rows, end
                if text.startswith(",[", pos + 1):
                    pos += 2
                    token_start = pos + 1
            else:
                token_start = pos + 1
        pos += 1
    raise ValueError("unterminated Maple Matrix")


def load_maple_lrp(path: str | Path, shape: tuple[int, int, int]) -> Scheme:
    path = Path(path)
    data = bz2.open(path, "rt").read() if path.suffix == ".bz2" else path.read_text()
    starts = [match.start() for match in _MATRIX_HEADER.finditer(data)]
    if len(starts) != 3:
        raise ValueError(f"expected exactly three LRP matrices, found {len(starts)}")
    parsed = [_parse_maple_matrix(data, start) for start in starts]
    m, n, p = shape
    rank = parsed[0][0]
    dimensions = tuple((item[0], item[1]) for item in parsed)
    expected = ((rank, m * n), (rank, n * p), (m * p, rank))
    if dimensions != expected:
        raise ValueError(f"LRP dimensions {dimensions}, expected {expected}")
    maple_u, maple_v = parsed[0][2], parsed[1][2]
    # Maple's LinearAlgebra vectorization is column-major for input matrices.
    U = [{(index % m) * n + index // m: value for index, value in row.items()} for row in maple_u]
    V = [{(index % n) * p + index // n: value for index, value in row.items()} for row in maple_v]
    # L and R are product-by-coordinate; P is coordinate-by-product.
    W: list[SparseRow] = [{} for _ in range(rank)]
    for coordinate, row in enumerate(parsed[2][2]):
        for product, value in row.items():
            W[product][coordinate] = value
    scheme = Scheme(
        field="Q", shape=shape, rank=rank, U=U, V=V, W=W, source=str(path),
        provenance={"format": "FMM-Lille Maple LRP", "sha256": sha256_file(path)},
        exact_coefficient_type="rational",
    )
    scheme.validate_dimensions()
    return scheme


def _flatten_matrix(rows: list[SparseRow], width: int) -> SparseRow:
    return {i * width + j: value for i, row in enumerate(rows) for j, value in row.items()}


def load_maple_tensor(path: str | Path, shape: tuple[int, int, int]) -> Scheme:
    """Load FMM-Lille's authoritative list of explicit ``Triad`` matrices."""
    path = Path(path)
    data = bz2.open(path, "rt").read() if path.suffix == ".bz2" else path.read_text()
    marker = data.find("Tensor:=TriadSet([")
    if marker < 0:
        raise ValueError("Tensor assignment is missing")
    body = data[marker:]
    starts = [match.start() for match in _MATRIX_HEADER.finditer(body)]
    if len(starts) % 3:
        raise ValueError(f"tensor contains {len(starts)} matrices, not complete triads")
    m, n, p = shape
    U: list[SparseRow] = []
    V: list[SparseRow] = []
    W: list[SparseRow] = []
    for product in range(len(starts) // 3):
        a = _parse_maple_matrix(body, starts[3 * product])
        b = _parse_maple_matrix(body, starts[3 * product + 1])
        c = _parse_maple_matrix(body, starts[3 * product + 2])
        if (a[0], a[1], b[0], b[1], c[0], c[1]) != (m, n, n, p, p, m):
            raise ValueError(f"triad {product} has incompatible matrix dimensions")
        U.append(_flatten_matrix(a[2], n))
        V.append(_flatten_matrix(b[2], p))
        # The third factor is explicitly a p×m matrix, i.e. C transpose.
        W.append({i * p + j: value for j, row in enumerate(c[2]) for i, value in row.items()})
    scheme = Scheme(
        field="Q", shape=shape, rank=len(U), U=U, V=V, W=W, source=str(path),
        provenance={"format": "FMM-Lille Maple Tensor/Triad", "sha256": sha256_file(path)},
        exact_coefficient_type="rational",
    )
    scheme.validate_dimensions()
    return scheme


def load_scheme(path: str | Path, shape: tuple[int, int, int] | None = None) -> Scheme:
    path = Path(path)
    if path.name.endswith((".mpl", ".mpl.bz2")):
        if shape is None:
            match = re.search(r"(\d+)x(\d+)x(\d+)", path.name)
            if not match:
                raise ValueError("shape is required for a Maple LRP file")
            shape = tuple(map(int, match.groups()))
        return load_maple_tensor(path, shape) if "_tensor.mpl" in path.name else load_maple_lrp(path, shape)
    return load_json_scheme(path)


def save_scheme(scheme: Scheme, path: str | Path) -> None:
    Path(path).write_text(json.dumps(scheme.to_json(), indent=2, sort_keys=True) + "\n")
