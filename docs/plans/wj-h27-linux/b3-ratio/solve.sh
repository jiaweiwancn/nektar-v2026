#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
stage=$1
name=$2
seed="$(realpath "$3")"
mode=$4
entry=$5
cooperative=$6
iterations=$7
radius=$8
frame="$(realpath "$9")"
shift 9
optimizer_seconds=$1
outer=$((cooperative+180))
case "$mode" in
    retained) inverse_estimate=0 ;;
    scalar-blocks) inverse_estimate=1750 ;;
    *) exit 2 ;;
esac
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
old="$(dirname "$r")/wj-c1-b3-bounded.xjpRylk7"
wrapper="$(dirname "$r")/wj-c1-steady-recovery.LOp8mfob/run-capped.sh"
case_dir="$r/$name"
test -s "$r/bin/$stage/finish.epoch"
bash capacity.sh "$name-workers.ps" 2
python3 "$suite/c1_steady_case.py" "$old/bounded-active/run.xml" "$case_dir" \
    --seed "$seed" --eta 0 --quadrature 10 16 --fd .000512 --factor-fd 4e-6 \
    --seconds "$cooperative" --newton "$iterations" --krylov 24 --cycles 1 \
    --inverse ColoredSparseV1 --inverse-mass .009539982703780335 --action-gap-audit --fixed-k-tangent \
    --steady-forcing .3 .001 .3 --merit-mass 0 --merit-reanchor 0 --merit-entry "$entry" --merit-work 96 \
    --maximum-inverse-age 20 --reuse-linear-iterations 24 --accepted-checkpoint "$outer" 110 "$inverse_estimate"
python3 - "$case_dir" "$suite/c1_bounded_active_optimizer.py" "$frame" "$radius" "$mode" "$optimizer_seconds" <<'R11_SESSION'
import json
from pathlib import Path
import sys

from c1_bounded_trials import load, physical_xml
from resources import identity
from task2r_verify import edit

case, operator, frame_path = map(Path,sys.argv[1:4])
radius, mode, seconds = sys.argv[4:]
session = case/"run.xml"
physical = physical_xml(session)
frame = load(frame_path)
assert frame["schema"]=="WJH27RecycledBasisV1" and 1<=frame["columns"]<=22
assert identity(Path(frame["input"]["path"]))==frame["input"]
script = identity(operator)
edit(session,{"WJDiagnosticMeritPolicy":"WJH27BoundedDirectionV2",
              "WJDiagnosticBoundedOptimizer":operator,"WJDiagnosticBoundedOptimizerSHA256":script["sha256"],
              "WJDiagnosticBoundedBasis":frame["input"]["path"],"WJDiagnosticBoundedBasisSHA256":frame["input"]["sha256"],
              "WJDiagnosticBoundedBasisBinding":frame["binding"],
              "WJDiagnosticBoundedResidualBlocks":"True" if mode=="scalar-blocks" else "False"},
             {"WJDiagnosticBoundedInitialRadius":radius,"WJDiagnosticBoundedComparisonRadius":.025,
              "WJDiagnosticBoundedMinimumRadius":1e-7,"WJDiagnosticBoundedMaximumRadius":.1,
              "WJDiagnosticBoundedColumns":frame["columns"],"WJDiagnosticBoundedProposals":4,
              "WJDiagnosticBoundedOptimizerIterations":1000,"WJDiagnosticBoundedOptimizerSeconds":seconds,
              "WJDiagnosticBoundedOptimizerRelativeGap":.1})
assert physical_xml(session)==physical
inputs=load(case/"inputs.json")
inputs[1]=identity(session)
inputs.extend((script,identity(frame_path),frame["input"]))
(case/"inputs.json").write_text(json.dumps(inputs,indent=2)+"\n")
R11_SESSION
date +%s > "$case_dir/start.epoch"
status=0
bash "$wrapper" "$case_dir" "$r/bin/$stage/IncNavierStokesSolver" "$outer" || status=$?
date +%s > "$case_dir/finish.epoch"
exit "$status"
