#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
export PYTHONPATH="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
while ! test -f bounded-main/exit-status || ! test -f bounded-active/exit-status; do
  test "$(($(date +%s)-$(cat start.epoch)))" -lt 8400
  sleep 12
done
remaining=$((10800-$(date +%s)+$(cat start.epoch)))
test "$remaining" -gt 1960
bash workers.sh final-audit-workers.ps 65
readarray -t selected < <(python3 - "$r" <<'R10_SELECT_ENDPOINT'
import json,sys
from pathlib import Path
from c1_bounded_trials import packed,save,load
from resources import identity
root=Path(sys.argv[1])
choices=[]
for stage in ("bounded-main","bounded-recycled","bounded-enriched","enriched-resume","bounded-active"):
    for path in (root/stage).glob("state.accepted-*.json"):
        p=load(path)
        assert p["regularization"]["eta"]==0 and p["iterate_custody"]=="LAST_ACCEPTED"
        choices.append((max(p["original_residuals"]),path))
norm, seed=min(choices,key=lambda x:x[0])
best=load(seed)
packed(best,best["state"]["state_sha256"])
mode="qualify" if norm<=1e-10 else "pending-all"
save(root/"selected-endpoint.json",{"checkpoint":identity(seed),"native_max":norm,
    "coefficient_sha256":best["state"]["state_sha256"],"mode":mode})
print(seed)
print(mode)
R10_SELECT_ENDPOINT
)
test "${#selected[@]}" = 2
bash run-endpoint.sh recycled "${selected[0]}" endpoint-rich "${selected[1]}" 1370
if test "$(cat endpoint-rich/admission.txt)" = QUALIFIED; then
  remaining=$((10800-$(date +%s)+$(cat start.epoch)))
  if test "$remaining" -gt 775; then
    bash run-endpoint.sh recycled "${selected[0]}" endpoint-rich paired 615
  fi
fi
date +%s > endpoint-rich/finish.epoch
