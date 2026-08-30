from germsynth.gf2_flip import (
    search_rank7,
    target_matrix_multiplication_tensor,
    tensor_of_scheme,
)
from germsynth.karatsuba_search import search_rank3
from germsynth.ternary_lift import lift_support_to_integers
from germsynth.tensor import (
    degree1_polynomial_multiplication_specification,
    matrix_multiplication_2x2_specification,
    sum_terms,
)


def test_karatsuba_is_discovered_exactly() -> None:
    result = search_rank3()
    specification = degree1_polynomial_multiplication_specification()
    assert len(result.terms) == 3
    assert sum_terms(result.terms, (2, 2, 3)) == specification.tensor


def test_flip_search_and_integer_lift() -> None:
    support = search_rank7(seed=8)
    assert len(support.scheme) == 7
    assert tensor_of_scheme(support.scheme) == target_matrix_multiplication_tensor()
    lifted = lift_support_to_integers(support.scheme)
    specification = matrix_multiplication_2x2_specification()
    assert len(lifted.terms) == 7
    assert sum_terms(lifted.terms, (4, 4, 4)) == specification.tensor
