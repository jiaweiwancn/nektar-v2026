#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
old="$(dirname "$r")/wj-c1-b3-bounded.xjpRylk7"
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
basis="$old/bounded-active/state.action-bounded-current-basis-base-3-proposal--1-sequence-9.json"
incumbent="$old/bounded-active/state.action-bounded-direction-base-3-proposal-1-sequence-11.json"
mkdir reduced
bash capacity.sh reduced/workers.ps 2
legacy="$old/optimizer-source/solvers/IncNavierStokesSolver/Tests/RANS/h27/c1_bounded_active_optimizer.py"
/usr/bin/time -v -o reduced/legacy-025.cost timeout 40s \
    python3 "$legacy" "$basis" .025 "$r/reduced/legacy-025"
for radius in .00625 .025 .003125; do
    /usr/bin/time -v -o "reduced/ratio-$radius.cost" timeout 40s \
        python3 "$suite/c1_bounded_active_optimizer.py" "$basis" "$radius" "$r/reduced/ratio-$radius" \
            --iterations 1000 --seconds 25 --relative-gap .1 --incumbent "$incumbent"
done
date +%s > reduced/finish.epoch
