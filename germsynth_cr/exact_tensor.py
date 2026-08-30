from __future__ import annotations

from germsynth_fusion.exact_tensor import audit_products, tensor_hash, verify_scheme as _verify
from germsynth_fusion.scheme_io import Scheme, SparseRow


SCREEN_PRIMES = (1000003, 1000033, 1000037)


def verify_scheme(scheme: Scheme):
    """Deterministic all-coordinate verification with three modular screens and exact Q/Z finish."""
    return _verify(scheme, SCREEN_PRIMES)


__all__ = ["Scheme", "SparseRow", "verify_scheme", "audit_products", "tensor_hash", "SCREEN_PRIMES"]
