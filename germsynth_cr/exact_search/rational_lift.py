from __future__ import annotations

from fractions import Fraction


def lift_coefficients(coefficients, modulus: int, bound: int = 1000):
    result = []
    for value in coefficients:
        centered = value if value <= modulus // 2 else value - modulus
        if abs(centered) > bound:
            return {"status": "FAIL", "reason": "centered residue exceeds lift bound"}
        result.append(Fraction(centered))
    return {"status": "CANDIDATE", "coefficients": result,
            "warning": "A lifted candidate is not valid until exact tensor verification passes."}
