from fractions import Fraction
from pathlib import Path

from germsynth_fusion.closure_extraction import component_residual, extract_components
from germsynth_fusion.exact_tensor import verify_scheme
from germsynth_fusion.scheme_io import load_scheme
from germsynth_fusion.scheme_io import Scheme
from germsynth_fusion.spectral_grammar import Edge, solve_critical_exponent
from decimal import Decimal


ROOT = Path(__file__).resolve().parents[1]
FMM = ROOT / "external" / "fmm-lille"


def scheme(shape):
    return load_scheme(FMM / shape / f"{shape}_tensor.mpl.bz2")


def test_small_exact_sources_pass():
    assert verify_scheme(scheme("2x3x3"))["status"] == "PASS"
    assert verify_scheme(scheme("4x9x10"))["status"] == "PASS"


def test_parent_source_has_direct_exact_counterexample():
    parent = scheme("8x27x30")
    a, b, c = 24, 729, 9
    actual = sum((u.get(a, Fraction(0)) * v.get(b, Fraction(0)) * w.get(c, Fraction(0))
                  for u, v, w in zip(parent.U, parent.V, parent.W)), Fraction(0))
    assert actual == 0  # target coefficient is 1
    check = verify_scheme(parent, (1000003,))
    assert check["status"] == "FAIL"
    assert check["modular_screens"][0]["residual_count"] == 28098


def test_automatic_census_recovers_nominal_14_but_no_exact_kernel():
    inner, outer, parent = scheme("2x3x3"), scheme("4x9x10"), scheme("8x27x30")
    extraction = extract_components(parent, outer, inner)
    reports = [component_residual(parent, outer, inner, component) for component in extraction.components]
    assert len(extraction.ordinary_slots) == 238
    assert len(extraction.ordinary_product_ids) == 3570
    assert sorted(report["product_count"] for report in reports) == [21, 29, 29, 29, 29, 29]
    assert sum(report["product_count"] for report in reports) + 3570 == 3736
    assert sum(report["fusion_gain"] for report in reports) == 14
    assert not any(report["standalone"] for report in reports)


def test_finite_field_final_verification_stays_in_field():
    # In GF(2), two identical scalar products cancel; over Q they would sum to 2.
    scheme = Scheme("GF(2)", (1, 1, 1), 3, [{0: Fraction(1)}] * 3,
                    [{0: Fraction(1)}] * 3, [{0: Fraction(1)}] * 3, "synthetic")
    assert verify_scheme(scheme)["status"] == "PASS"
    scheme.field = "Q"
    assert verify_scheme(scheme)["status"] == "FAIL"


def test_periodic_two_type_spectral_radius():
    edges = [Edge(2, 0, 1, Decimal(2)), Edge(8, 1, 0, Decimal(2))]
    result = solve_critical_exponent(edges, 2)
    assert abs(float(result["gamma"]) - 2.0) < 1e-10
