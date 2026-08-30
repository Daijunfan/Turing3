from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .exact_tensor import verify_scheme
from .scheme_io import Scheme, SparseRow


def _kron_row(outer: SparseRow, inner: SparseRow, outer_shape: tuple[int, int], inner_shape: tuple[int, int]) -> SparseRow:
    orows, ocols = outer_shape
    irows, icols = inner_shape
    width = ocols * icols
    result: SparseRow = {}
    for oi_flat, ov in outer.items():
        oi, oj = divmod(oi_flat, ocols)
        for ii_flat, iv in inner.items():
            ii, ij = divmod(ii_flat, icols)
            result[(oi * irows + ii) * width + (oj * icols + ij)] = ov * iv
    return result


def kronecker_scheme(outer: Scheme, inner: Scheme) -> Scheme:
    om, on, op = outer.shape
    im, inn, ip = inner.shape
    U, V, W = [], [], []
    for slot in range(outer.rank):
        for product in range(inner.rank):
            U.append(_kron_row(outer.U[slot], inner.U[product], (om, on), (im, inn)))
            V.append(_kron_row(outer.V[slot], inner.V[product], (on, op), (inn, ip)))
            W.append(_kron_row(outer.W[slot], inner.W[product], (om, op), (im, ip)))
    return Scheme(
        field=outer.field, shape=(om * im, on * inn, op * ip), rank=len(U), U=U, V=V, W=W,
        source=f"Kronecker({outer.source},{inner.source})",
        provenance={"constructor": "plain Kronecker", "outer_rank": outer.rank, "inner_rank": inner.rank},
        exact_coefficient_type="rational",
    )


@dataclass
class FusionConstructor:
    organ_certificate: dict

    def match(self, slot_shapes, field: str, orientation=(0, 1, 2)) -> bool:
        return (all(tuple(shape) == tuple(self.organ_certificate["slot_shape"]) for shape in slot_shapes)
                and field == self.organ_certificate["field"] and sorted(orientation) == [0, 1, 2]
                and self.organ_certificate.get("status") == "PASS")

    def instantiate(self, slot_embeddings):
        if self.organ_certificate.get("status") != "PASS":
            raise ValueError("cannot instantiate an unverified fusion organ")
        return {"certificate": self.organ_certificate, "slot_embeddings": slot_embeddings}

    def replace_independent_slots(self, scheme: Scheme, *_args, **_kwargs) -> Scheme:
        if self.organ_certificate.get("status") != "PASS":
            raise ValueError("replacement refused: fusion certificate is not PASS")
        raise NotImplementedError("verified coordinate rebasing is required for this organ")

    def emit_exact_scheme(self, scheme: Scheme) -> Scheme:
        return scheme

    def verify(self, scheme: Scheme):
        return verify_scheme(scheme)
