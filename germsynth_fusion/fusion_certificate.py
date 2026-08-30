from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from .closure_extraction import component_residual
from .exact_tensor import tensor_hash, verify_scheme
from .fusion_hypergraph import HypergraphComponent
from .scheme_io import Scheme, SparseRow


def encode_row(row: SparseRow) -> list[list[int]]:
    return [[index, value.numerator, value.denominator] for index, value in sorted(row.items())]


def encode_products(scheme: Scheme, product_ids: list[int] | None = None) -> list[dict[str, Any]]:
    ids = list(range(scheme.rank)) if product_ids is None else product_ids
    return [{"id": product, "U": encode_row(scheme.U[product]), "V": encode_row(scheme.V[product]),
             "W": encode_row(scheme.W[product])} for product in ids]


def certificate_hash(certificate: dict[str, Any]) -> str:
    payload = dict(certificate)
    payload.pop("certificate_sha256", None)
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def parent_certificate(scheme: Scheme, verification: dict[str, Any], extraction_log: dict[str, Any]) -> dict[str, Any]:
    certificate = {
        "format": "germsynth-fusion-parent-certificate-v1",
        "status": verification["status"],
        "field": scheme.field,
        "shape": list(scheme.shape),
        "rank": scheme.rank,
        "exact_coefficient_type": scheme.exact_coefficient_type,
        "vectorization": scheme.vectorization,
        "products": encode_products(scheme),
        "exact_tensor_hash": tensor_hash(scheme),
        "source": scheme.source,
        "provenance": scheme.provenance,
        "verification": verification,
        "automatic_extraction_log": extraction_log,
    }
    certificate["certificate_sha256"] = certificate_hash(certificate)
    return certificate


def organ_certificate(
    parent: Scheme, outer: Scheme, inner: Scheme, component: HypergraphComponent, kron_ranks
) -> dict[str, Any]:
    residual = component_residual(parent, outer, inner, component)
    certificate = {
        "format": "germsynth-fusion-organ-certificate-v1",
        "status": "PASS" if residual["standalone"] else "FAIL",
        "field": parent.field,
        "ambient_shape": list(parent.shape),
        "outer_shape": list(outer.shape),
        "slot_shape": list(inner.shape),
        "slot_count": len(component.slot_ids),
        "rank": len(component.product_ids),
        "products": encode_products(parent, component.product_ids),
        "slot_ids": component.slot_ids,
        "slot_embeddings": [
            {"slot_id": slot, "U": encode_row(outer.U[slot]), "V": encode_row(outer.V[slot]),
             "W": encode_row(outer.W[slot])} for slot in component.slot_ids
        ],
        "product_slot_incidence": [],
        "kron_ranks": [list(kron_ranks[product]) for product in component.product_ids],
        "mixed_product_ids": [product for product in component.product_ids
                              if kron_ranks[product] != (1, 1, 1)],
        "plain_product_ids": [product for product in component.product_ids
                              if kron_ranks[product] == (1, 1, 1)],
        "fusion_gain": residual["fusion_gain"],
        "standalone": residual["standalone"],
        "dependency_closure": residual["classification"],
        "exact_tensor_certificate": residual,
        "parent_scheme_provenance": parent.provenance,
    }
    certificate["certificate_sha256"] = certificate_hash(certificate)
    return certificate


def write_certificate(certificate: dict[str, Any], path: str | Path) -> None:
    Path(path).write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
