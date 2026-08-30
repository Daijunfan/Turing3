"""Scale-weighted spectral certificates for multi-type algorithm germs.

For recursive calls e: i -> j that shrink linear problem scale by s_e > 1,
define A_ij(gamma) = sum_e multiplicity_e * s_e**(-gamma).  The critical
exponent is the positive root rho(A(gamma)) = 1.  A one-type, uniform germ
reduces to gamma = log_b(rank).
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np


@dataclass(frozen=True, slots=True)
class ScaleTransition:
    source: int
    target: int
    shrink_factor: float
    multiplicity: int = 1

    def __post_init__(self) -> None:
        if self.source < 0 or self.target < 0:
            raise ValueError("type indices must be nonnegative")
        if self.shrink_factor <= 1:
            raise ValueError("shrink_factor must exceed one")
        if self.multiplicity <= 0:
            raise ValueError("multiplicity must be positive")


@dataclass(frozen=True, slots=True)
class ScaleWeightedGerm:
    type_count: int
    transitions: tuple[ScaleTransition, ...]

    def __post_init__(self) -> None:
        if self.type_count <= 0:
            raise ValueError("type_count must be positive")
        if any(
            transition.source >= self.type_count or transition.target >= self.type_count
            for transition in self.transitions
        ):
            raise ValueError("transition type index is out of range")

    def operator(self, exponent: float) -> np.ndarray:
        matrix = np.zeros((self.type_count, self.type_count), dtype=float)
        for transition in self.transitions:
            matrix[transition.source, transition.target] += (
                transition.multiplicity * transition.shrink_factor ** (-exponent)
            )
        return matrix

    def spectral_radius(self, exponent: float) -> float:
        eigenvalues = np.linalg.eigvals(self.operator(exponent))
        return float(max(abs(value) for value in eigenvalues))

    def critical_exponent(self, tolerance: float = 1e-12) -> float:
        if not self.transitions:
            return 0.0
        low, high = 0.0, 1.0
        if self.spectral_radius(low) < 1.0:
            return 0.0
        while self.spectral_radius(high) > 1.0:
            high *= 2.0
            if high > 1e6 or not isfinite(high):
                raise RuntimeError("failed to bracket critical exponent")
        for _ in range(200):
            middle = (low + high) / 2.0
            if self.spectral_radius(middle) > 1.0:
                low = middle
            else:
                high = middle
            if high - low <= tolerance:
                break
        return (low + high) / 2.0
