from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class ContactClass(str, Enum):
    INDEPENDENT = "INDEPENDENT"
    KRONECKER = "KRONECKER"
    AXIS_CONCAT = "AXIS_CONCAT"
    SERENDIPITOUS = "SERENDIPITOUS"
    PAN_PAIR = "PAN_PAIR"
    DISJOINT_SUM = "DISJOINT_SUM"
    GAUGE_EQUIVALENT = "GAUGE_EQUIVALENT"
    GENUINE_MIXED_CONTACT = "GENUINE_MIXED_CONTACT"
    UNCLASSIFIED = "UNCLASSIFIED"


@dataclass
class ContactOrgan:
    slot_ids: tuple[int, ...]
    rank: int
    baseline_cost: int
    classification: ContactClass
    exact: bool
    provenance: dict = field(default_factory=dict)

    @property
    def contact_gain(self) -> int:
        return self.baseline_cost - self.rank

    @property
    def is_edge(self) -> bool:
        return self.exact and self.contact_gain > 0


@dataclass
class ContactCertificate:
    field: str
    commutative: bool
    organ: ContactOrgan
    tensor_hash: str


VERIFIED_EXPLICIT = "VERIFIED_EXPLICIT"
VERIFIED_DERIVED = "VERIFIED_DERIVED"
METADATA_ONLY = "METADATA_ONLY"
INVALID = "INVALID"
REPAIRED = "REPAIRED"
QUARANTINED = "QUARANTINED"
