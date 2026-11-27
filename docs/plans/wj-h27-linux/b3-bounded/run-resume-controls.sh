#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
checkpoint="$r/bounded-main/state.accepted-000001.json"
while ! test -s "$checkpoint"; do
  if test -s "$r/bounded-main/exit-status"; then exit 4; fi
  test "$(($(date +%s)-$(cat start.epoch)))" -lt 10200
  sleep 12
done
bash workers.sh resume-controls-workers.ps 66
# A short cooperative limit triggers the inverse-reserve guard after restoring
# the controller and before preparing a new cold inverse.
for entry in Resume Reanchor; do
  radius=.025
  if [ "$entry" = Reanchor ]; then radius=.0125; fi
  label="cold-${entry,,}"
  bash run-bounded.sh bounded "$label" "$checkpoint" "$entry" 64 1960 1 "$radius" 16 22 \
    > "$label.launch.log" 2>&1 || test -s "$label/state.diagnostic-steady.json"
done
if bash run-bounded.sh bounded changed-resume "$checkpoint" Resume 64 1960 1 .0125 16 22 \
    > changed-resume.launch.log 2>&1; then exit 7; fi
grep -m1 "original-merit resume requires the exact accepted coefficients, policy and equation/source bindings" changed-resume/run.log
python3 - "$checkpoint" <<'VERIFY_RESUME_R10'
import hashlib,json,sys
from pathlib import Path
q=Path(sys.argv[1])
a, same, changed=(json.loads(p.read_text()) for p in
    (q,Path("cold-resume/state.diagnostic-steady.json"),Path("cold-reanchor/state.diagnostic-steady.json")))
assert a["state"]["state_sha256"]==same["state"]["state_sha256"]==changed["state"]["state_sha256"]
assert a["original_residuals"]==same["original_residuals"]==changed["original_residuals"]
old=a["merit_controller"]
new=same["merit_controller"]
fresh=changed["merit_controller"]
assert new["state"]==old["state"] and new["policy_sha256"]==old["policy_sha256"]
assert new["accepted_displacement_sha256"]==old["accepted_displacement_sha256"]
assert new["accepted_displacement"]==old["accepted_displacement"]
assert new["entry"]=="Resume" and fresh["entry"]=="Reanchor" and fresh["policy_sha256"]!=old["policy_sha256"]
assert fresh["state"]["radius"]["next"]==.0125 and not fresh.get("accepted_displacement")
for case in ("cold-resume","cold-reanchor"):
    ledger=json.loads((Path(case)/"state.steady.json").read_text())
    assert not ledger["accepted_newton_steps"] and not ledger["physical_time_advanced"]
    assert not list(Path(case).glob("state.accepted-*.json"))
    assert "before inverse" in ledger["reason"]
assert not (Path("changed-resume")/"state.diagnostic-steady.json").exists()
def pin(p):
    return {"path":str(p.resolve()),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"bytes":p.stat().st_size}
with Path("resume-controls.json").open("x") as stream:
    json.dump({"seed":pin(q),"unchanged_policy":pin(Path("cold-resume/state.diagnostic-steady.json")),
               "explicit_reanchor":pin(Path("cold-reanchor/state.diagnostic-steady.json")),
               "mutation_rejection":pin(Path("changed-resume/run.log")),
               "restored_policy_sha256":new["policy_sha256"],"state":new["state"],
               "accepted_displacement_sha256":new["accepted_displacement_sha256"],
               "reanchored_state":fresh["state"]},stream,indent=2)
    stream.write("\n")
VERIFY_RESUME_R10
