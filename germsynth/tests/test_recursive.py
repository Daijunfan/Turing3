import random

from germsynth.germ import BilinearGerm
from germsynth.gf2_flip import search_rank7
from germsynth.karatsuba_search import search_rank3
from germsynth.recursive import (
    OperationCounter,
    germ_matrix_multiply,
    germ_polynomial_multiply,
    naive_matrix_multiply,
    naive_polynomial_multiply,
)
from germsynth.tensor import (
    degree1_polynomial_multiplication_specification,
    matrix_multiplication_2x2_specification,
)
from germsynth.ternary_lift import lift_support_to_integers


def _germs() -> tuple[BilinearGerm, BilinearGerm]:
    karatsuba = search_rank3()
    matrix_support = search_rank7(seed=8)
    matrix_lift = lift_support_to_integers(matrix_support.scheme)
    return (
        BilinearGerm(
            "Karatsuba germ",
            degree1_polynomial_multiplication_specification(),
            2,
            karatsuba.terms,
        ),
        BilinearGerm(
            "rank-7 matrix germ",
            matrix_multiplication_2x2_specification(),
            2,
            matrix_lift.terms,
        ),
    )


def test_recursive_matrix_algorithm() -> None:
    _, germ = _germs()
    rng = random.Random(11)
    for n in (1, 2, 4, 8, 16):
        left = [[rng.randint(-3, 3) for _ in range(n)] for _ in range(n)]
        right = [[rng.randint(-3, 3) for _ in range(n)] for _ in range(n)]
        counter = OperationCounter()
        assert germ_matrix_multiply(left, right, germ, counter) == naive_matrix_multiply(left, right)
        assert counter.scalar_multiplications == 7 ** (n.bit_length() - 1)


def test_recursive_polynomial_algorithm() -> None:
    germ, _ = _germs()
    rng = random.Random(12)
    for n in (1, 2, 4, 8, 16, 32, 64):
        left = [rng.randint(-5, 5) for _ in range(n)]
        right = [rng.randint(-5, 5) for _ in range(n)]
        counter = OperationCounter()
        assert germ_polynomial_multiply(left, right, germ, counter) == naive_polynomial_multiply(left, right)
        assert counter.scalar_multiplications == 3 ** (n.bit_length() - 1)
