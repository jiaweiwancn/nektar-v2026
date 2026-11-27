#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
stage=$1
case_name=$2
seed=$3
entry=$4
cooperative=$5
outer=$6
iterations=$7
initial_radius=$8
columns=$9
inverse_age=${10}
frame_manifest=${11}
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
operator="$r/optimizer-source/solvers/IncNavierStokesSolver/Tests/RANS/h27/c1_bounded_active_optimizer.py"
export PYTHONPATH="$suite"
test -s "$r/bin/$stage/finish.epoch"
bash workers.sh "$r/$case_name-workers.ps" 67
old="/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/hybrid-fine/run.xml"
case_dir="$r/$case_name"
python3 "$suite/c1_steady_case.py" "$old" "$case_dir" --seed "$seed" --eta 0 --quadrature 10 16 \
  --fd .000512 --factor-fd 4e-6 --seconds "$cooperative" --newton "$iterations" --krylov "$columns" --cycles 1 \
  --inverse ColoredSparseV1 --inverse-mass .009539982703780335 --action-gap-audit --fixed-k-tangent \
  --steady-forcing .3 .001 .3 --merit-mass 0 --merit-reanchor 0 --merit-entry "$entry" --merit-work 96 \
  --maximum-inverse-age "$inverse_age" --reuse-linear-iterations "$columns" --accepted-checkpoint "$outer" 110 0
python3 - "$case_dir" "$operator" "$initial_radius" "$columns" "$frame_manifest" <<'RADIUS_OPTIONS'
from pathlib import Path
from task2r_verify import edit
from resources import identity
import json,sys
case,operator=map(Path,sys.argv[1:3])
session=case/"run.xml"
script=identity(operator)
frame=json.loads(Path(sys.argv[5]).read_text())
assert frame["schema"]=="WJH27RecycledBasisV1" and frame["columns"]==int(sys.argv[4])
assert identity(Path(frame["input"]["path"]))==frame["input"]
edit(session,{"WJDiagnosticMeritPolicy":"WJH27BoundedDirectionV1","WJDiagnosticBoundedOptimizer":operator,
              "WJDiagnosticBoundedOptimizerSHA256":script["sha256"],"WJDiagnosticBoundedBasis":frame["input"]["path"],
              "WJDiagnosticBoundedBasisSHA256":frame["input"]["sha256"],
              "WJDiagnosticBoundedBasisBinding":frame["binding"]},
             {"WJDiagnosticBoundedInitialRadius":sys.argv[3],"WJDiagnosticBoundedMinimumRadius":1e-7,
              "WJDiagnosticBoundedMaximumRadius":.1,"WJDiagnosticBoundedColumns":sys.argv[4],
              "WJDiagnosticBoundedProposals":4})
path=case/"inputs.json"
inputs=json.loads(path.read_text())
inputs[1]=identity(session)
inputs.extend((script, identity(Path(sys.argv[5])), frame["input"]))
path.write_text(json.dumps(inputs,indent=2)+"\n")
RADIUS_OPTIONS
date +%s > "$case_dir/start.epoch"
bash "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" \
  "$case_dir" "$r/bin/$stage/IncNavierStokesSolver" "$outer"
date +%s > "$case_dir/finish.epoch"
