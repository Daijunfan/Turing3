#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p build results
exec > >(tee results/full_reproduction.log) 2>&1

PY="$ROOT/germsynth/.venv/bin/python"
if [[ ! -x "$PY" ]]; then
  python3 -m venv "$ROOT/germsynth/.venv"
  "$PY" -m pip install -q -e "$ROOT/germsynth" pytest
fi

echo '[1/6] Existing GermSynth-R/F reproduction'
FUSION_LOG=results/fusion_reproduction.log ./reproduce_fusion.sh

echo '[2/6] Source forensics, independent parsers, residual autopsy and contact compiler'
"$PY" run_contact_repair.py

echo '[3/6] Python tests'
"$PY" -m pytest -q tests
(cd germsynth && "$PY" -m pytest -q)

echo '[4/6] Compile independent C++17 verifiers'
for source in verify_8x27x30_3744 verify_residual_completion verify_contact_certificate verify_novel_candidate; do
  c++ -std=c++17 -O2 -Wall -Wextra -pedantic "independent_verifiers/${source}.cpp" -o "build/${source}"
done

echo '[5/6] Independent C++ controls and required rejections'
./build/verify_contact_certificate certificates/known_4x4x4_rank48.json
if ./build/verify_8x27x30_3744 certificates/matmul_8x27x30_rank3744.json; then
  echo 'unverified rank-3744 certificate unexpectedly passed' >&2; exit 1
fi
if ./build/verify_residual_completion certificates/residual_completion.json; then
  echo 'missing residual completion unexpectedly passed' >&2; exit 1
fi
if ./build/verify_novel_candidate certificates/novel_candidates/manifest.json; then
  echo 'empty novel manifest unexpectedly passed' >&2; exit 1
fi

echo '[6/6] Final machine statuses'
jq -r '.statuses | to_entries[] | "\(.key)=\(.value)"' results/latest_summary.json
echo 'FULL_CONTACT_REPAIR_REPRODUCTION=FAIL'
echo 'REPRODUCTION_PIPELINE_COMPLETED=TRUE'
