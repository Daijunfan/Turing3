import json

from germsynth.certificates import write_certificate
from germsynth.germ import BilinearGerm
from germsynth.karatsuba_search import search_rank3
from germsynth.tensor import degree1_polynomial_multiplication_specification


def test_certificate_is_self_consistent(tmp_path) -> None:
    result = search_rank3()
    germ = BilinearGerm(
        "Karatsuba germ",
        degree1_polynomial_multiplication_specification(),
        2,
        result.terms,
    )
    path = write_certificate(germ, tmp_path / "certificate.json")
    payload = json.loads(path.read_text())
    assert payload["local_identity_verified"] is True
    assert payload["rank"] == 3
    assert payload["target_tensor_sha256"] == payload["reconstructed_tensor_sha256"]
