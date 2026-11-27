#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
export PYTHONPATH="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
bash workers.sh finite-negative-workers.ps 66
python3 - "$r" <<'NEGATIVE_XML'
from copy import deepcopy
import hashlib,json,sys
from pathlib import Path
from shutil import copyfile
from task2r_verify import edit
root=Path(sys.argv[1])
original=root/"frozen-trials"
metadata=json.loads((original/"frozen-trial-inputs.json").read_text())
for name in ("stale-base","stale-scale","stale-vector"):
    case=root/name
    case.mkdir()
    session=case/"run.xml"
    copyfile(original/"run.xml",session)
    options={"WJOutputPrefix":case/"state"}
    if name=="stale-base":
        seed=Path("/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/curvature/state.accepted-000001.json")
        data=seed.read_bytes()
        options.update(WJInitialFields=seed,WJInitialSourceIdentity=json.loads(data)["state"]["source_identity"],
                       WJDiagnosticSeedSHA256=hashlib.sha256(data).hexdigest())
    elif name=="stale-scale":
        info=metadata["fixed_base"]
        binding=hashlib.sha256("|".join(["WJH27FrozenFiniteBindingV1",info["base"],info["ordering"],
                          "0"*64,info["row_scale"],metadata["input"]["sha256"]]).encode()).hexdigest()
        options["WJDiagnosticActionProbeBinding"]=binding
    else:
        payload=bytearray(Path(metadata["input"]["path"]).read_bytes())
        payload[0]^=1
        bad=case/"changed-direction.bin"
        bad.write_bytes(payload)
        options["WJDiagnosticActionProbe"]=bad
    edit(session,options,{})
NEGATIVE_XML
pids=()
for case in stale-base stale-scale stale-vector; do
  (if bash "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" \
      "$r/$case" "$r/bin/finite/IncNavierStokesSolver" 190; then
      exit 90
   fi
   if [ "$case" = stale-vector ]; then
      grep -m1 'fixed direction size or digest differs' "$r/$case/run.log"
   else
      grep -m1 'finite-trial input base, free ordering or scales differ' "$r/$case/run.log"
   fi
   test -z "$(find "$r/$case" -maxdepth 1 -name "state.action-*.json" -print -quit)"
   ) &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
