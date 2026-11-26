#!/usr/bin/env bash
set -euo pipefail
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
old="/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/hybrid-fine"
basis=("$r"/basis/state.action-frozen-size-basis-*.json)
test "${#basis[@]}" -eq 1
test -s "${basis[0]}"
bash workers.sh margin-workers.ps 68
mkdir margin
/usr/bin/time -v -o margin/optimizer.cost python3 frozen_size.py propose "${basis[0]}" \
    "$old/state.action-merit-unresolved-base-1-proposal-0-sequence-3.json" "$r/margin" \
    --radii 12000 12500 13000
bash verify.sh "$r/margin"
