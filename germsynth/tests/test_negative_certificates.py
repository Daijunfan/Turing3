import json
from pathlib import Path
import subprocess
import sys


def test_independent_verifier_rejects_corrupted_germ(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    payload = json.loads((root / "certificates" / "matrix_rank7_germ.json").read_text())
    payload["terms"][0]["w"][0] += 1
    corrupted = tmp_path / "corrupted_germ.json"
    corrupted.write_text(json.dumps(payload))
    completed = subprocess.run(
        [sys.executable, str(root / "independent_verify.py"), str(corrupted)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode != 0


def test_independent_pool_verifier_rejects_corrupted_phenotype(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    payload = json.loads((root / "certificates" / "matrix_phenotype_pool.json").read_text())
    payload["phenotypes"][0]["terms"][0]["u"][0] += 1
    corrupted = tmp_path / "corrupted_pool.json"
    corrupted.write_text(json.dumps(payload))
    completed = subprocess.run(
        [sys.executable, str(root / "independent_verify_pool.py"), str(corrupted)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode != 0
