r"""Exact bilinear tensor primitives used by GermSynth.

No floating point arithmetic is used in synthesis or verification.  A bilinear
algorithm is represented by a sum of rank-one integer tensors

    T = sum_s u_s \otimes v_s \otimes w_s.

The coordinates of ``u`` and ``v`` select linear forms in the two inputs and
``w`` specifies how the resulting scalar product contributes to outputs.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import product
from typing import Iterable, Sequence

Vector = tuple[int, ...]
Tensor = tuple[int, ...]


@dataclass(frozen=True, slots=True)
class RankOneTerm:
    u: Vector
    v: Vector
    w: Vector

    def tensor(self) -> Tensor:
        return outer(self.u, self.v, self.w)

    def to_dict(self) -> dict[str, list[int]]:
        return {"u": list(self.u), "v": list(self.v), "w": list(self.w)}

    @staticmethod
    def from_dict(data: dict[str, Sequence[int]]) -> "RankOneTerm":
        return RankOneTerm(tuple(data["u"]), tuple(data["v"]), tuple(data["w"]))


@dataclass(frozen=True, slots=True)
class BilinearSpecification:
    name: str
    a_dim: int
    b_dim: int
    out_dim: int
    tensor: Tensor
    kind: str
    base_shape: tuple[int, ...]

    def __post_init__(self) -> None:
        expected = self.a_dim * self.b_dim * self.out_dim
        if len(self.tensor) != expected:
            raise ValueError(f"tensor has {len(self.tensor)} entries, expected {expected}")

    def tensor_hash(self) -> str:
        return tensor_hash(self.tensor)


def outer(u: Sequence[int], v: Sequence[int], w: Sequence[int]) -> Tensor:
    return tuple(a * b * c for a in u for b in v for c in w)


def tensor_add(a: Sequence[int], b: Sequence[int]) -> Tensor:
    if len(a) != len(b):
        raise ValueError("tensor length mismatch")
    return tuple(x + y for x, y in zip(a, b, strict=True))


def tensor_sub(a: Sequence[int], b: Sequence[int]) -> Tensor:
    if len(a) != len(b):
        raise ValueError("tensor length mismatch")
    return tuple(x - y for x, y in zip(a, b, strict=True))


def sum_terms(terms: Iterable[RankOneTerm], dims: tuple[int, int, int]) -> Tensor:
    size = dims[0] * dims[1] * dims[2]
    result = [0] * size
    for term in terms:
        if (len(term.u), len(term.v), len(term.w)) != dims:
            raise ValueError("rank-one term dimension mismatch")
        for idx, value in enumerate(term.tensor()):
            result[idx] += value
    return tuple(result)


def tensor_hash(tensor: Sequence[int]) -> str:
    payload = ",".join(str(x) for x in tensor).encode("ascii")
    return sha256(payload).hexdigest()


def matrix_multiplication_2x2_specification() -> BilinearSpecification:
    """The exact 4 x 4 x 4 tensor for multiplying two 2x2 matrices.

    Input block order is row-major: 00, 01, 10, 11.  Output order is the same.
    """
    a_dim = b_dim = out_dim = 4
    target = [0] * (a_dim * b_dim * out_dim)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ai = i * 2 + k
                bj = k * 2 + j
                ci = i * 2 + j
                target[(ai * b_dim + bj) * out_dim + ci] += 1
    return BilinearSpecification(
        name="2x2 matrix multiplication",
        a_dim=a_dim,
        b_dim=b_dim,
        out_dim=out_dim,
        tensor=tuple(target),
        kind="matrix_multiplication",
        base_shape=(2, 2, 2),
    )


def degree1_polynomial_multiplication_specification() -> BilinearSpecification:
    """Tensor for (a0+a1*x)(b0+b1*x) -> (c0,c1,c2)."""
    a_dim, b_dim, out_dim = 2, 2, 3
    target = [0] * (a_dim * b_dim * out_dim)
    for i, j, k in ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 2)):
        target[(i * b_dim + j) * out_dim + k] += 1
    return BilinearSpecification(
        name="degree-1 polynomial multiplication",
        a_dim=a_dim,
        b_dim=b_dim,
        out_dim=out_dim,
        tensor=tuple(target),
        kind="polynomial_multiplication",
        base_shape=(2, 2, 3),
    )


def evaluate_bilinear_terms(
    terms: Sequence[RankOneTerm], a: Sequence[int], b: Sequence[int]
) -> tuple[int, ...]:
    if not terms:
        return ()
    out_dim = len(terms[0].w)
    out = [0] * out_dim
    for term in terms:
        if len(term.u) != len(a) or len(term.v) != len(b) or len(term.w) != out_dim:
            raise ValueError("evaluation dimension mismatch")
        left = sum(c * x for c, x in zip(term.u, a, strict=True))
        right = sum(c * x for c, x in zip(term.v, b, strict=True))
        value = left * right
        for idx, coefficient in enumerate(term.w):
            out[idx] += coefficient * value
    return tuple(out)


def evaluate_specification(
    specification: BilinearSpecification, a: Sequence[int], b: Sequence[int]
) -> tuple[int, ...]:
    if len(a) != specification.a_dim or len(b) != specification.b_dim:
        raise ValueError("specification evaluation dimension mismatch")
    out = [0] * specification.out_dim
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            for k in range(specification.out_dim):
                coefficient = specification.tensor[
                    (i * specification.b_dim + j) * specification.out_dim + k
                ]
                out[k] += coefficient * ai * bj
    return tuple(out)


def exhaustive_value_vectors(dim: int, alphabet: Sequence[int]) -> Iterable[tuple[int, ...]]:
    return product(alphabet, repeat=dim)
