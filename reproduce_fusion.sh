#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p build results
exec > >(tee results/full_reproduction.log) 2>&1

BASE="$ROOT/germsynth"
PY="$BASE/.venv/bin/python"
if [[ ! -x "$PY" ]]; then
  python3 -m venv "$BASE/.venv"
  "$BASE/.venv/bin/python" -m pip install -q -e "$BASE" pytest
fi

echo '[1/7] GermSynth-R baseline regeneration'
(
  cd "$BASE"
  export PYTHONPATH="$BASE"
  "$PY" run_all.py
  "$PY" -m pytest -q
  "$PY" independent_verify.py certificates/karatsuba_germ.json
  "$PY" independent_verify.py certificates/matrix_rank7_germ.json
  "$PY" independent_verify_pool.py certificates/matrix_phenotype_pool.json
  "$PY" emit_cpp.py certificates/matrix_rank7_germ.json generated/matrix_rank7_germ.cpp
  c++ -std=c++17 -O2 -Wall -Wextra -pedantic generated/matrix_rank7_germ.cpp -o build/matrix_rank7_germ_verify
  ./build/matrix_rank7_germ_verify
  "$PY" emit_cpp_pool.py certificates/matrix_phenotype_pool.json generated/matrix_pool_verify.cpp
  c++ -std=c++17 -O2 -Wall -Wextra -pedantic generated/matrix_pool_verify.cpp -o build/matrix_pool_verify
  ./build/matrix_pool_verify
)

echo '[2/7] Exact Python extraction and certificates'
"$PY" run_fusion.py

echo '[3/7] GermSynth-F Python tests'
"$PY" -m pytest -q tests

echo '[4/7] Compile independent C++17 verifiers'
c++ -std=c++17 -O2 -Wall -Wextra -pedantic independent_verifiers/verify_parent_scheme.cpp -o build/verify_parent_scheme
c++ -std=c++17 -O2 -Wall -Wextra -pedantic independent_verifiers/verify_fusion_kernel.cpp -o build/verify_fusion_kernel

echo '[5/7] Positive controls and coefficient-mutation rejection'
./build/verify_parent_scheme certificates/reference_2x3x3_rank15.json
if ./build/verify_parent_scheme certificates/reference_2x3x3_rank15.json --mutate-first; then
  echo 'mutation unexpectedly passed parent verifier' >&2
  exit 1
fi
./build/verify_fusion_kernel certificates/reference_plain_233_rank15.json
if ./build/verify_fusion_kernel certificates/reference_plain_233_rank15.json --mutate-first; then
  echo 'mutation unexpectedly passed kernel verifier' >&2
  exit 1
fi

echo '[6/7] Required invalid-source certificates must be rejected'
if ./build/verify_parent_scheme certificates/fusion_parent_8x27x30_rank3736.json; then
  echo 'invalid 3736 parent unexpectedly passed' >&2
  exit 1
fi
if ./build/verify_fusion_kernel certificates/fusion_2x233_rank29.json; then
  echo 'invalid rank-29 kernel unexpectedly passed' >&2
  exit 1
fi

echo '[7/7] Final status and package'
echo 'BASELINE_REPRODUCTION=PASS'
echo 'PARENT_EXACT_VERIFICATION=FAIL'
echo 'AUTOMATIC_SLOT_INFERENCE=PASS'
echo 'FUSION_KERNEL_EXTRACTION=FAIL'
echo 'FULL_GAIN_14_EXPLAINED=FAIL'
echo 'CXX_INDEPENDENT_VERIFICATION=FAIL'
echo 'CONSTRUCTOR_PASS=FAIL'
echo 'ALL_SIX_DEFICITS_CLOSED=FAIL'
echo 'NOVEL_FUSION_KERNEL=FAIL'
echo 'NEW_EXACT_RANK=FAIL'
echo 'NEW_EXPONENT=FAIL'
echo 'TURING_PATH_PASS=FAIL'
echo 'REPRODUCTION_COMPLETED=TRUE (required invalid claims were rejected)'

zip -qr GermSynth-F_result.zip \
  germsynth_fusion independent_verifiers tests results certificates external/SOURCES.json external/fmm-lille \
  README_FUSION.md reproduce_fusion.sh run_fusion.py \
  -x '*/__pycache__/*'
