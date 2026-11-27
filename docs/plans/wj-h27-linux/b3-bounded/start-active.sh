#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
bash workers.sh active-model-workers.ps 66
git clone --shared --no-checkout "$r/source" optimizer-source > optimizer-sync.log 2>&1
git -C optimizer-source fetch "$r/active-helper.bundle" HEAD >> optimizer-sync.log 2>&1
git -C optimizer-source checkout --detach FETCH_HEAD >> optimizer-sync.log 2>&1
git -C optimizer-source rev-parse HEAD > active-helper-source.txt
sha256sum "$r/active-helper.bundle" > active-helper-bundle.sha256
test -z "$(git -C optimizer-source status --porcelain)"
optimizer="$r/optimizer-source/solvers/IncNavierStokesSolver/Tests/RANS/h27/c1_bounded_active_optimizer.py"
basis=("$r"/enriched-resume/state.action-bounded-current-basis-base-0-proposal--1-sequence-1.json)
test -s "${basis[0]}"
/usr/bin/time -v -o active-model.cost timeout --kill-after=2s 120s python3 "$optimizer" "${basis[0]}" .025 active-fixed-model
python3 - <<'ACTIVATION_R10'
from pathlib import Path
import json
d=json.loads(Path("active-fixed-model.json").read_text())
print(d)
assert d["size"]>0 and d["model_score"]<1., "active-gradient proposal did not find an improving model"
ACTIVATION_R10
seed="$r/bounded-enriched/state.accepted-000002.json"
frame="$r/recycled-frames/enriched.json"
columns=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["columns"])' "$frame")
bash run-active.sh recycled bounded-active "$seed" Reanchor 2350 2560 12 .025 "$columns" 22 "$frame"
