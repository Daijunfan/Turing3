from __future__ import annotations

import bz2
import re
from fractions import Fraction
from pathlib import Path

from .exact_tensor import Scheme, SparseRow


INPUT_TERM = re.compile(r"([+-]?)(?:(\d+(?:/\d+)?)\*)?([AB])_(\d+)_(\d+)")
OUTPUT_TERM = re.compile(r"([+-]?)(?:(\d+(?:/\d+)?)\*)?m_(\d+)")


def _linear_form(expression: str, symbol: str, rows: int, columns: int) -> SparseRow:
    matches = list(INPUT_TERM.finditer(expression))
    if INPUT_TERM.sub("", expression).replace("(", "").replace(")", ""):
        raise ValueError(f"unsupported raw input expression: {expression}")
    result: SparseRow = {}
    for match in matches:
        if match.group(3) != symbol:
            raise ValueError("raw expression mixes input matrices")
        sign = -1 if match.group(1) == "-" else 1
        value = sign * Fraction(match.group(2) or 1)
        i, j = int(match.group(4)) - 1, int(match.group(5)) - 1
        if not (0 <= i < rows and 0 <= j < columns):
            raise ValueError("raw input coordinate out of range")
        coordinate = i * columns + j
        result[coordinate] = result.get(coordinate, Fraction()) + value
    return {coordinate: value for coordinate, value in result.items() if value}


def parse_raw(path: str | Path) -> Scheme:
    """Symbolically execute FMM-Lille MUL/ADD without using Tensor/LRP parsing code."""
    path = Path(path)
    text = bz2.open(path, "rt").read() if path.suffix == ".bz2" else path.read_text()
    if "MUL:=[" not in text:
        return _parse_procedure(path, text)
    declarations = [tuple(map(int, pair)) for pair in
                    re.findall(r"^[ABC]:=Matrix\((\d+),\s*(\d+),", text, re.MULTILINE)]
    if len(declarations) != 3:
        raise ValueError("raw A/B/C declarations missing")
    (m, n), (n2, p), (m2, p2) = declarations
    if n != n2 or (m, p) != (m2, p2):
        raise ValueError("raw matrix dimensions are incompatible")
    mul_start = text.index("MUL:=[") + len("MUL:=[")
    mul_end = text.index("]:\nADD:=", mul_start)
    multiplication_lines = [line.strip().rstrip(",") for line in text[mul_start:mul_end].splitlines() if line.strip()]
    U: list[SparseRow] = []
    V: list[SparseRow] = []
    for expected, line in enumerate(multiplication_lines, 1):
        match = re.fullmatch(r"m_(\d+)=\((.*)\)\*\((.*)\)", line)
        if not match or int(match.group(1)) != expected:
            raise ValueError(f"raw multiplication sequence breaks at {expected}")
        U.append(_linear_form(match.group(2), "A", m, n))
        V.append(_linear_form(match.group(3), "B", n, p))
    W: list[SparseRow] = [{} for _ in U]
    add_start = text.index("ADD:=[") + len("ADD:=[")
    add_end = text.index("]:\nmap(expand", add_start)
    outputs_seen = set()
    for line in text[add_start:add_end].splitlines():
        line = line.strip().rstrip(",")
        if not line or line == "0=0":
            continue
        equation = re.fullmatch(r"C_(\d+)_(\d+)=(.*)", line)
        if not equation:
            raise ValueError(f"unsupported raw output equation: {line[:80]}")
        i, j = int(equation.group(1)) - 1, int(equation.group(2)) - 1
        output = i * p + j
        outputs_seen.add(output)
        expression = equation.group(3)
        terms = list(OUTPUT_TERM.finditer(expression))
        if OUTPUT_TERM.sub("", expression):
            raise ValueError(f"unsupported raw output expression C_{i+1}_{j+1}")
        for term in terms:
            sign = -1 if term.group(1) == "-" else 1
            value = sign * Fraction(term.group(2) or 1)
            product = int(term.group(3)) - 1
            if not 0 <= product < len(W):
                raise ValueError("raw output references undefined multiplication")
            W[product][output] = W[product].get(output, Fraction()) + value
    if outputs_seen != set(range(m * p)):
        raise ValueError("raw program does not define every output")
    scheme = Scheme("Q", (m, n, p), len(U), U, V, W, str(path),
                    {"parser": "germsynth_cr.maple_raw_parser", "representation": "symbolic MUL/ADD"})
    scheme.validate_dimensions()
    return scheme


class _Value:
    def __init__(self, kind: str, data=None):
        self.kind = kind
        self.data = data if data is not None else {}


def _combine(left: _Value, right: _Value, sign=1) -> _Value:
    if left.kind == "S" and left.data == 0:
        return _scale(right, sign)
    if right.kind == "S" and right.data == 0:
        return left
    if left.kind != right.kind:
        raise ValueError(f"raw SLP adds incompatible {left.kind}/{right.kind}")
    if left.kind == "S":
        return _Value("S", left.data + sign * right.data)
    data = dict(left.data)
    for key, value in right.data.items():
        updated = data.get(key, Fraction()) + sign * value
        if updated:
            data[key] = updated
        else:
            data.pop(key, None)
    return _Value(left.kind, data)


def _scale(value: _Value, scalar) -> _Value:
    scalar = Fraction(scalar)
    if value.kind == "S":
        return _Value("S", value.data * scalar)
    return _Value(value.kind, {key: coefficient * scalar for key, coefficient in value.data.items() if coefficient * scalar})


TOKENS = re.compile(r"[A-Za-z]\w*(?:\[\d+,\d+\])?|\d+(?:/\d+)?|[-+*()]")


class _Expression:
    def __init__(self, expression: str, environment: dict[str, _Value], products):
        self.tokens = TOKENS.findall(expression.replace(" ", ""))
        if "".join(self.tokens) != expression.replace(" ", ""):
            raise ValueError(f"unsupported raw SLP expression {expression}")
        self.position = 0
        self.environment = environment
        self.products = products

    def parse(self):
        value = self._sum()
        if self.position != len(self.tokens):
            raise ValueError("raw SLP expression has trailing tokens")
        return value

    def _sum(self):
        value = self._product()
        while self.position < len(self.tokens) and self.tokens[self.position] in ("+", "-"):
            operator = self.tokens[self.position]
            self.position += 1
            value = _combine(value, self._product(), -1 if operator == "-" else 1)
        return value

    def _product(self):
        value = self._atom()
        while self.position < len(self.tokens) and self.tokens[self.position] == "*":
            self.position += 1
            right = self._atom()
            if value.kind == "S":
                value = _scale(right, value.data)
            elif right.kind == "S":
                value = _scale(value, right.data)
            elif {value.kind, right.kind} == {"A", "B"}:
                a, b = (value, right) if value.kind == "A" else (right, value)
                product = len(self.products)
                self.products.append((a.data, b.data))
                value = _Value("P", {product: Fraction(1)})
            else:
                raise ValueError("raw SLP contains non-bilinear multiplication")
        return value

    def _atom(self):
        token = self.tokens[self.position]
        self.position += 1
        if token == "-":
            return _scale(self._atom(), -1)
        if token == "+":
            return self._atom()
        if token == "(":
            value = self._sum()
            if self.tokens[self.position] != ")":
                raise ValueError("raw SLP parenthesis mismatch")
            self.position += 1
            return value
        if token[0].isdigit():
            return _Value("S", Fraction(token))
        if token not in self.environment:
            raise ValueError(f"raw SLP references undefined {token}")
        return self.environment[token]


def _parse_procedure(path: Path, text: str) -> Scheme:
    shape_match = re.search(r"FastMatMul_(\d+)x(\d+)x(\d+):=proc", text)
    if not shape_match:
        raise ValueError("raw procedure shape missing")
    m, n, p = map(int, shape_match.groups())
    environment: dict[str, _Value] = {}
    for i in range(m):
        for j in range(n):
            environment[f"A[{i+1},{j+1}]"] = _Value("A", {i * n + j: Fraction(1)})
    for i in range(n):
        for j in range(p):
            environment[f"B[{i+1},{j+1}]"] = _Value("B", {i * p + j: Fraction(1)})
    products: list[tuple[SparseRow, SparseRow]] = []
    outputs: dict[int, _Value] = {}
    for statement in text.split(";"):
        statement = statement.strip()
        if ":=" not in statement or statement.startswith("FastMatMul_"):
            continue
        left, expression = map(str.strip, statement.split(":=", 1))
        if expression.startswith("Matrix("):
            continue
        call = re.fullmatch(r"LinearAlgebra:-MatrixMatrixMultiply\((\w+),(\w+)\)", expression)
        if call:
            a, b = environment[call.group(1)], environment[call.group(2)]
            if (a.kind, b.kind) != ("A", "B"):
                raise ValueError("raw MatrixMatrixMultiply operands are not A/B linear")
            product = len(products)
            products.append((a.data, b.data))
            value = _Value("P", {product: Fraction(1)})
        else:
            value = _Expression(expression, environment, products).parse()
        output_match = re.fullmatch(r"C\[(\d+),(\d+)\]", left)
        if output_match:
            i, j = int(output_match.group(1)) - 1, int(output_match.group(2)) - 1
            if value.kind != "P":
                raise ValueError("raw output is not bilinear")
            outputs[i * p + j] = value
        else:
            environment[left] = value
    if set(outputs) != set(range(m * p)):
        raise ValueError("raw SLP does not define every output")
    U = [pair[0] for pair in products]
    V = [pair[1] for pair in products]
    W: list[SparseRow] = [{} for _ in products]
    for output, value in outputs.items():
        for product, coefficient in value.data.items():
            W[product][output] = coefficient
    scheme = Scheme("Q", (m, n, p), len(products), U, V, W, str(path),
                    {"parser": "germsynth_cr.maple_raw_parser", "representation": "symbolic Maple SLP"})
    scheme.validate_dimensions()
    return scheme
