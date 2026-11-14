#!/usr/bin/env bash
# Linux only. Reuse the pinned runtime; every output and history is fresh.
# The physical pair is EXPECTED to exit 1 on temporal:Pi. Do not bypass it.
set -euo pipefail
archive='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-8-20261114'
test "$#" -eq 1 || { echo 'usage: reproduce-strip.sh NEW_ABSOLUTE_DIRECTORY' >&2; exit 2; }
out="$1"
case "$out" in /*) ;; *) echo 'output must be absolute' >&2; exit 2;; esac
test ! -e "$out" || { echo 'use a new output directory' >&2; exit 2; }
mkdir "$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$archive/tmp" OMPI_MCA_orte_tmpdir_base=/tmp/nektar-ct8-ipc
runtime="$archive/build-v3/runtime"
export LD_LIBRARY_PATH="$runtime/lib"
source_root="$archive/source"
tests="$source_root/solvers/IncNavierStokesSolver/Tests/RANS"
contract="$source_root/docs/plans/cao-tamura-accepted-8-evidence/contract.json"
data="$source_root/docs/plans/cao-tamura-validation-data"
guard() { local tag="$1"; shift; python3 -B "$tests/ct_accepted3_command.py" --control "$out/launch-$tag" -- "$@"; }
guard reference python3 -B "$tests/ct_accepted8_reference.py" --data "$data" --contract "$contract" --runtime "$runtime" --output "$out/reference"
guard preflight python3 -B "$tests/ct_accepted8_native.py" --mode preflight --runtime "$runtime" --reference "$out/reference/results.json" --contract "$contract" --output "$out/preflight"
guard compatibility python3 -B "$tests/ct_accepted8_controls.py" --runtime "$runtime" --reference "$out/reference/results.json" --output "$out/compatibility"
guard physical python3 -B "$tests/ct_accepted8_native.py" --mode advance --runtime "$runtime" --reference "$out/reference/results.json" --contract "$contract" --preflight "$out/preflight" --controls "$out/compatibility" --target .008 --output "$out/T008"
