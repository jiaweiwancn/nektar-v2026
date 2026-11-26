#!/usr/bin/env bash
set -euo pipefail
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
old="/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/hybrid-fine"
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
if [ "$#" -gt 0 ]; then
  case_dir="$(realpath "$1")"
else
  case_dir="$r"
fi
seed="$old/state.accepted-000001.json"
bash workers.sh "$case_dir/verify-workers.ps" 68
test -s "$case_dir/candidates.json"
python3 "$suite/c1_steady_case.py" "$old/run.xml" "$case_dir/verify" --seed "$seed" --eta 0 --quadrature 10 16 \
 --evaluate-only --seconds 1050 --fd .000512 --steady-forcing .3 .001 .3
python3 - "$case_dir" <<'PY'
from pathlib import Path
from task2r_verify import edit
from resources import identity
import json,sys
root=Path(sys.argv[1])
case=root/"verify"
session=case/"run.xml"
manifest=json.loads((root/"candidates.json").read_text())
matrix=Path(manifest["input"]["path"])
assert identity(matrix)["sha256"]==manifest["input"]["sha256"]
edit(session, {"WJDiagnosticActionProbe":matrix,
               "WJDiagnosticActionProbeSHA256":manifest["input"]["sha256"], "WJDiagnosticFixedKTangent":"True"},
              {"WJDiagnosticActionProbeCount":len(manifest["candidates"])})
path=case/"inputs.json"
inputs=json.loads(path.read_text())
inputs[1]=identity(session)
inputs.append(identity(matrix))
path.write_text(json.dumps(inputs,indent=2)+"\n")
PY
date +%s > "$case_dir/verify-start.epoch"
bash "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" "$case_dir/verify" "$r/build/IncNavierStokesSolver" 1170
/usr/bin/time -v -o "$case_dir/summarize.cost" python3 frozen_size.py summarize "$case_dir" "$case_dir/verify"
