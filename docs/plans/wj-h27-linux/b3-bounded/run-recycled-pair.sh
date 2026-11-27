#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
bash workers.sh recycled-pair-workers.ps 63
seed="$(find bounded-main -maxdepth 1 -name 'state.accepted-*.json' -print | sort | tail -1)"
test -s "$seed"
readarray -t settings < <(python3 - "$seed" <<'PAIR_SEED'
from pathlib import Path
import json,sys
path=Path(sys.argv[1])
x=json.loads(path.read_text())
assert max(x["original_residuals"])>1e-10
assert x["merit_controller"]["state"]["accepted"]>=3
print(x["accepted_checkpoint"]["sequence"]-1)
print(x["merit_controller"]["state"]["radius"]["next"])
PAIR_SEED
)
basis=("$r"/bounded-main/state.action-bounded-current-basis-base-"${settings[0]}"-*.json)
test "${#basis[@]}" = 1 && test -s "${basis[0]}"
python3 "$suite/c1_recycled_basis.py" "$seed" "${basis[0]}" "$r/recycled-frames"
python3 - "$seed" "${basis[0]}" "${settings[1]}" <<'PAIR_SETUP'
import json, sys
from pathlib import Path
from resources import identity
p=Path("recycled-comparison-inputs.json")
with p.open("x") as s:
    json.dump({"seed":identity(sys.argv[1]),"basis":identity(sys.argv[2]),
               "initial_radius":float(sys.argv[3]),"frame_control":identity("recycled-frames/recycled.json"),
               "frame_enriched":identity("recycled-frames/enriched.json"),
               "maximum_native_seconds":640,"accepted_iteration_bound":2},s,indent=2)
    s.write("\n")
PAIR_SETUP
pids=()
for frame in recycled enriched; do
  manifest="$r/recycled-frames/$frame.json"
  cols="$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["columns"])' "$manifest")"
  (bash run-recycled.sh recycled "bounded-$frame" "$r/$seed" Reanchor 640 810 2 "${settings[1]}" "$cols" 22 "$manifest" \
    > "bounded-$frame.launch.log" 2>&1 || test -s "bounded-$frame/state.diagnostic-steady.json") &
  pids+=("$!")
done
statuses=0
for child in "${pids[@]}"; do wait "$child" || statuses=1; done
exit "$statuses"
