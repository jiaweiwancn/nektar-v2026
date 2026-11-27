#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
stage=$1
name=$2
shift 2
old="$(dirname "$r")/wj-c1-b3-bounded.xjpRylk7"
export PYTHONPATH="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
wrapper="$(dirname "$r")/wj-c1-steady-recovery.LOp8mfob/run-capped.sh"
exe="$r/bin/$stage/IncNavierStokesSolver"
test -s "$r/bin/$stage/finish.epoch"
bash capacity.sh "$name-workers.ps" 2
case_dir="$r/$name"
manifest="$r/source/docs/plans/wj-h27-c1-b3-bounded-nonlinear-progress-evidence.json"
/usr/bin/time -v -o "$name.setup.cost" python3 replay_ratio.py "$case_dir" \
    --old "$old" --manifest "$manifest" "$@"
bash "$wrapper" "$case_dir" "$exe" 660 || test "$(cat "$case_dir/exit-status")" = 2
python3 replay_ratio.py "$case_dir"
date +%s > "$case_dir/finish.epoch"
