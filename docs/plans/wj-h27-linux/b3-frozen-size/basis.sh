#!/usr/bin/env bash
set -euo pipefail
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
old="/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/hybrid-fine"
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
seed="$old/state.accepted-000001.json"
printf '%s  %s\n' 84dae45ef585f24321137adb4b131ffc3157c15644ed1cd16ac7221b34caabb6 "$seed" d75919798347090e68e573f59b5c5be99d1446bd392feefc2ca278db1a2506da "$old/run.xml" | sha256sum -c -
python3 "$suite/c1_steady_case.py" "$old/run.xml" "$r/basis" \
 --seed "$seed" --eta 0 --quadrature 10 16 --fd .000512 --factor-fd 4e-6 \
 --seconds 2150 --newton 1 --krylov 16 --cycles 1 \
 --inverse ColoredSparseV1 --inverse-mass .009539982703780335 --action-gap-audit --fixed-k-tangent \
 --steady-forcing .3 .001 .3 --merit-mass 0 --merit-reanchor 0 --merit-entry Reanchor --merit-work 96
python3 - "$r/basis" <<'PY'
from pathlib import Path
from task2r_verify import edit
from resources import identity
import json,sys
case=Path(sys.argv[1])
session=case/"run.xml"
edit(session, {}, {"WJDiagnosticFrozenSizeColumns":16})
path=case/"inputs.json"
inputs=json.loads(path.read_text())
table=Path(inputs[2]["path"])
if identity(table)["sha256"]!="3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1":
    raise SystemExit("reference table digest differs")
inputs[1]=identity(session)
path.write_text(json.dumps(inputs,indent=2)+"\n")
PY
bash workers.sh basis-workers.ps 68
date +%s > basis-start.epoch
bash "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" "$r/basis" "$r/build/IncNavierStokesSolver" 2360
