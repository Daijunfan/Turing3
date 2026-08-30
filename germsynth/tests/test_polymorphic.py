import random

from germsynth.germ import BilinearGerm
from germsynth.gf2_flip import search_rank7
from germsynth.polymorphic import PhenotypePool, polymorphic_matrix_multiply, random_selector
from germsynth.recursive import OperationCounter, naive_matrix_multiply
from germsynth.tensor import matrix_multiplication_2x2_specification
from germsynth.ternary_lift import lift_support_to_integers


def test_arbitrary_nodewise_phenotype_mixture() -> None:
    specification = matrix_multiplication_2x2_specification()
    germs = []
    for seed in range(4):
        support = search_rank7(seed)
        lift = lift_support_to_integers(support.scheme)
        germs.append(BilinearGerm(f"phenotype-{seed}", specification, 2, lift.terms))
    pool = PhenotypePool(tuple(germs))
    rng = random.Random(13)
    n = 16
    left = [[rng.randint(-2, 2) for _ in range(n)] for _ in range(n)]
    right = [[rng.randint(-2, 2) for _ in range(n)] for _ in range(n)]
    counter = OperationCounter()
    actual = polymorphic_matrix_multiply(
        left, right, pool, random_selector(len(germs), 99), counter
    )
    assert actual == naive_matrix_multiply(left, right)
    assert counter.scalar_multiplications == 7 ** 4


def test_exact_regenerative_cover() -> None:
    from germsynth.regenerative_cover import RegenerativeCoverOptimizer

    specification = matrix_multiplication_2x2_specification()
    unique = {}
    for seed in range(20):
        support = search_rank7(seed)
        lift = lift_support_to_integers(support.scheme)
        germ = BilinearGerm(f"phenotype-{seed}", specification, 2, lift.terms, provenance={"seed": seed})
        key = tuple((term.u, term.v, term.w) for term in germ.terms)
        unique.setdefault(key, germ)
    pool = PhenotypePool(tuple(unique.values()))
    optimizer = RegenerativeCoverOptimizer(pool)
    cover = optimizer.minimum_single_fault_cover_with_tiebreak()
    assert len(cover.phenotype_indices) == 3
    assert cover.covered_failure_sets == cover.total_failure_sets
