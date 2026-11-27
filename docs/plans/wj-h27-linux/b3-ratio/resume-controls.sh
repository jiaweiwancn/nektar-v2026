#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
seed="$r/retained/state.accepted-000001.json"
test -s "$seed"
regular=0
if test -s resume-control/exit-status; then
    regular=$(cat resume-control/exit-status)
else
    bash solve.sh ratio resume-control "$seed" retained Resume 18 1 .00625 "$r/frame/recycled.json" 12 ||
        regular=$?
fi
case "$regular" in
    1|2) test -s resume-control/state.steady.json ;;
    *) exit "$regular" ;;
esac
changed=0
if test -s changed-policy/exit-status; then
    changed=$(cat changed-policy/exit-status)
else
    bash solve.sh ratio changed-policy "$seed" retained Resume 24 1 .003125 "$r/frame/recycled.json" 12 ||
        changed=$?
fi
test "$changed" != 0
rg_bin=$(command -v rg || true)
if test -n "$rg_bin"; then
    "$rg_bin" -q 'original-merit resume requires the exact accepted coefficients, policy and equation/source bindings' changed-policy/run.log
else
    grep -q 'original-merit resume requires the exact accepted coefficients, policy and equation/source bindings' changed-policy/run.log
fi
python3 - "$seed" <<'R11_CHECK_RESUME'
import sys
from pathlib import Path
from c1_bounded_trials import load, packed, save
from resources import identity

seed=Path(sys.argv[1])
before=load(seed)
after=load(Path("resume-control/state.diagnostic-steady.json"))
ledger=load(Path("resume-control/state.steady.json"))
assert ledger["accepted_newton_steps"]==0 and ledger["reason"].startswith("DEADLINE:")
assert not list(Path("resume-control").glob("state.accepted-*.json"))
assert not list(Path("changed-policy").glob("state.accepted-*.json"))
packed(after,before["state"]["state_sha256"])
assert before["original_residuals"]==after["original_residuals"]
a,b=(p["merit_controller"] for p in (before,after))
assert a["entry"]=="Reanchor" and b["entry"]=="Resume"
assert {k:v for k,v in a.items() if k!="entry"}=={k:v for k,v in b.items() if k!="entry"}
save(Path("resume-controls.json"),{"start":identity(seed),
    "restored":identity(Path("resume-control/state.diagnostic-steady.json")),
    "controller":a["state"],"accepted_displacement_sha256":a["accepted_displacement_sha256"],
    "policy_sha256":a["policy_sha256"],"changed_log":identity(Path("changed-policy/run.log")),
    "recorded_entry_transition":[a["entry"],b["entry"]],"unchanged_custody":True})
R11_CHECK_RESUME
