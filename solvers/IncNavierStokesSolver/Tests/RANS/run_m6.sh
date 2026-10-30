#!/usr/bin/env bash
# Reproduce the M6 evidence sequentially. Study jobs may use up to 16 workers;
# every MPI rank and numerical library uses one computational thread.
set -euo pipefail
if [ "$(uname -s)" != Linux ] || [ "$#" -lt 2 ]; then
  printf 'Usage on Linux: run_m6.sh BUILD_DIRECTORY RUN_DIRECTORY [STUDY_JOBS=2]\n' >&2
  exit 2
fi
RANS_BUILD=$(realpath "$1")
RANS_RUN=$(realpath -m "$2")
RANS_JOBS=${3:-2}
case "$RANS_JOBS" in 1|2|3|4) ;; *) printf 'Use 1-4 study jobs (at most 16 ranks).\n' >&2; exit 2;; esac
RANS_TESTS=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RANS_RUN/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 TMPDIR="$RANS_RUN/tmp"
RANS_PHASE=''
trap 'code=$?; if [ -n "$RANS_PHASE" ]; then printf "%s\n" "$code" > "$RANS_RUN/final-$RANS_PHASE.status"; fi' EXIT
finish_phase() {
  printf '0\n' > "$RANS_RUN/final-$RANS_PHASE.status"
  printf 'Phase %s passed\n' "$RANS_PHASE"
  RANS_PHASE=''
}
RANS_PHASE=B
ctest --test-dir "$RANS_BUILD" --parallel 1 --output-on-failure > "$RANS_RUN/final-ctest.log" 2>&1
python3 "$RANS_TESTS/validate_m6.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-verification" --suite pilots --jobs 1 > "$RANS_RUN/final-pilots.log" 2>&1
python3 "$RANS_TESTS/validate_m6.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-verification" --suite restart --jobs 1 > "$RANS_RUN/final-restart.log" 2>&1
finish_phase
RANS_PHASE=E
python3 "$RANS_TESTS/verify_outflow_m6.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-outflow-equivalence" > "$RANS_RUN/final-outflow-equivalence.log" 2>&1
finish_phase
RANS_PHASE=F
python3 "$RANS_TESTS/validate_m6.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-standard" --suite standard --jobs 1 > "$RANS_RUN/final-standard.log" 2>&1
finish_phase
RANS_PHASE=A
python3 "$RANS_TESTS/validate_m6.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-study" --suite study --jobs "$RANS_JOBS" > "$RANS_RUN/final-study.log" 2>&1
finish_phase
RANS_PHASE=C
python3 "$RANS_TESTS/validate_m6.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-columns" --study-output "$RANS_RUN/final-study" --suite columns --jobs 1 > "$RANS_RUN/final-columns.log" 2>&1
python3 "$RANS_TESTS/validate_m6.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-robust" --suite robust --jobs 1 > "$RANS_RUN/final-robust.log" 2>&1
python3 "$RANS_TESTS/plot_m6.py" --study "$RANS_RUN/final-study" --output "$RANS_RUN/final-figures" > "$RANS_RUN/final-plots.log" 2>&1
finish_phase
RANS_PHASE=D
for milestone in 3 4 5; do
  python3 "$RANS_TESTS/validate_m${milestone}.py" --build "$RANS_BUILD" --output "$RANS_RUN/final-regressions/m${milestone}" > "$RANS_RUN/final-regression-m${milestone}.log" 2>&1
done
finish_phase
python3 "$RANS_TESTS/collect_m6.py" --run "$RANS_RUN" --build "$RANS_BUILD" --output "$RANS_RUN/final-figures/evidence.json"
