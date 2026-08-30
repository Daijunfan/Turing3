from math import isclose, log2, sqrt

from germsynth.spectral import ScaleTransition, ScaleWeightedGerm


def test_one_type_exponents() -> None:
    strassen = ScaleWeightedGerm(1, (ScaleTransition(0, 0, 2.0, 7),))
    karatsuba = ScaleWeightedGerm(1, (ScaleTransition(0, 0, 2.0, 3),))
    assert isclose(strassen.critical_exponent(), log2(7), abs_tol=1e-10)
    assert isclose(karatsuba.critical_exponent(), log2(3), abs_tol=1e-10)


def test_multitype_scale_operator() -> None:
    germ = ScaleWeightedGerm(
        2,
        (
            ScaleTransition(0, 0, 2.0, 2),
            ScaleTransition(0, 1, 2.0, 1),
            ScaleTransition(1, 0, 2.0, 1),
        ),
    )
    expected = log2(1 + sqrt(2))
    assert isclose(germ.critical_exponent(), expected, abs_tol=1e-10)
