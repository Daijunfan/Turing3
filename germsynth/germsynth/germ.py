"""Proof-carrying bilinear algorithm germs."""
from __future__ import annotations

from dataclasses import dataclass, field
from math import log
from typing import Any, Sequence

from .tensor import (
    BilinearSpecification,
    RankOneTerm,
    evaluate_bilinear_terms,
    evaluate_specification,
    sum_terms,
    tensor_hash,
)


@dataclass(frozen=True, slots=True)
class BilinearGerm:
    name: str
    specification: BilinearSpecification
    block_factor: int
    terms: tuple[RankOneTerm, ...]
    coefficient_domain: str = "integers"
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def rank(self) -> int:
        return len(self.terms)

    @property
    def complexity_exponent(self) -> float:
        if self.block_factor <= 1:
            raise ValueError("block_factor must be greater than one")
        return log(self.rank, self.block_factor)

    def reconstructed_tensor(self) -> tuple[int, ...]:
        spec = self.specification
        return sum_terms(self.terms, (spec.a_dim, spec.b_dim, spec.out_dim))

    def verify_local_identity(self) -> bool:
        return self.reconstructed_tensor() == self.specification.tensor

    def verify_exhaustively(self, alphabet: Sequence[int]) -> tuple[bool, int]:
        """Evaluate every pair of base inputs over ``alphabet`` exactly."""
        from itertools import product

        checked = 0
        for a in product(alphabet, repeat=self.specification.a_dim):
            for b in product(alphabet, repeat=self.specification.b_dim):
                checked += 1
                actual = evaluate_bilinear_terms(self.terms, a, b)
                expected = evaluate_specification(self.specification, a, b)
                if actual != expected:
                    return False, checked
        return True, checked

    def certificate_payload(self) -> dict[str, Any]:
        spec = self.specification
        reconstructed = self.reconstructed_tensor()
        return {
            "schema": "germsynth-bilinear-germ-v1",
            "name": self.name,
            "coefficient_domain": self.coefficient_domain,
            "problem": {
                "kind": spec.kind,
                "name": spec.name,
                "base_shape": list(spec.base_shape),
                "a_dim": spec.a_dim,
                "b_dim": spec.b_dim,
                "out_dim": spec.out_dim,
            },
            "block_factor": self.block_factor,
            "rank": self.rank,
            "terms": [term.to_dict() for term in self.terms],
            "target_tensor": list(spec.tensor),
            "target_tensor_sha256": spec.tensor_hash(),
            "reconstructed_tensor_sha256": tensor_hash(reconstructed),
            "local_identity_verified": reconstructed == spec.tensor,
            "complexity_certificate": {
                "recurrence": f"T(n) = {self.rank} T(n/{self.block_factor}) + lower-order work",
                "scalar_multiplication_exponent": self.complexity_exponent,
                "closed_form_on_powers": f"M({self.block_factor}^k) = {self.rank}^k",
            },
            "provenance": self.provenance,
        }
