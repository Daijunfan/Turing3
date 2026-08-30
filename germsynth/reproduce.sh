#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p build generated results

export PYTHONPATH="$ROOT"

printf '\n[1/8] Exact search, lifting, recursive and polymorphic experiments\n'
python3 run_all.py

printf '\n[2/8] Unit and negative-certificate tests\n'
python3 -m pytest -q

printf '\n[3/8] Independent Python germ verification\n'
python3 independent_verify.py certificates/karatsuba_germ.json
python3 independent_verify.py certificates/matrix_rank7_germ.json

printf '\n[4/8] Independent Python pool/cover verification\n'
python3 independent_verify_pool.py certificates/matrix_phenotype_pool.json \
  | tee results/pool_independent_verification.txt

CXX="${CXX:-g++}"

printf '\n[5/8] Emit and compile standalone C++ germ checker\n'
python3 emit_cpp.py certificates/matrix_rank7_germ.json generated/matrix_rank7_germ.cpp
"$CXX" -std=c++17 -O2 -Wall -Wextra -pedantic \
  generated/matrix_rank7_germ.cpp -o build/matrix_rank7_germ_verify
./build/matrix_rank7_germ_verify | tee results/cpp_verification.txt

printf '\n[6/8] Emit and compile standalone C++ phenotype-pool checker\n'
python3 emit_cpp_pool.py certificates/matrix_phenotype_pool.json generated/matrix_pool_verify.cpp
"$CXX" -std=c++17 -O2 -Wall -Wextra -pedantic \
  generated/matrix_pool_verify.cpp -o build/matrix_pool_verify
./build/matrix_pool_verify | tee results/cpp_pool_verification.txt

printf '\n[7/8] Build source manifest\n'
find . -type f \
  ! -path './.pytest_cache/*' \
  ! -path './*/__pycache__/*' \
  ! -path './build/*' \
  ! -name 'MANIFEST.sha256' \
  ! -name 'full_reproduction.log' \
  -print0 | sort -z | xargs -0 sha256sum > MANIFEST.sha256

printf '\n[8/8] Final status\n'
grep -q '"status": "PASS"' results/results.json
grep -q 'CXX_INDEPENDENT_VERIFICATION=PASS' results/cpp_verification.txt
grep -q 'CXX_POOL_INDEPENDENT_VERIFICATION=PASS' results/cpp_pool_verification.txt
printf 'FULL_REPRODUCTION=PASS\n'
