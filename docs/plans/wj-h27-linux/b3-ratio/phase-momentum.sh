#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
export PYTHONPATH="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
while ! test -s momentum-decision.json; do
    remaining=$((10800-$(date +%s)+$(cat start.epoch)))
    bash capacity.sh momentum-selection-workers.ps 1 || { sleep 15; continue; }
    python3 - "$remaining" <<'R11_MOMENTUM_PICK'
from pathlib import Path
import sys
from c1_bounded_trials import load, packed, save
from resources import identity

directory=Path("scalars")
checkpoints=sorted(directory.glob("state.accepted-*.json"))
ledger=directory/"state.steady.json"
terminal=load(ledger)["reason"] if ledger.exists() else None
remaining=int(sys.argv[1])
decision=None
report={"scalar_termination":terminal,"remaining_seconds":remaining}
if remaining<=1850+1820:
    decision="INSUFFICIENT_FACTOR_AND_AUDIT_BUDGET"
elif checkpoints:
    seed=checkpoints[-1]
    packet=load(seed)
    coefficients=packet["state"]["state_sha256"]
    packed(packet,coefficients)
    blocks=packet["original_residuals"]
    report.update(packet=identity(seed.resolve()),coefficients=coefficients,original_blocks=blocks)
    stalled=bool(terminal and terminal.startswith(("BOUNDED_STAGNATION:",
        "BOUNDED_RESIDUAL_SPAN_EXHAUSTED:","BOUNDED_INVERSE_AGE:")))
    slow_momentum=False
    if len(checkpoints)>=6 and max(blocks)==blocks[0]:
        samples=[load(path)["original_residuals"] for path in checkpoints[-3:]]
        relative=(samples[0][0]-samples[-1][0])/samples[0][0]
        slow_momentum=all(max(row)==row[0] for row in samples) and relative<.03
        report["last_three_accepted_original_blocks"]=samples
        report["last_two_relative_momentum_decrease"]=relative
    if max(blocks)<=1.e-10:
        decision="NATIVE_ROOT_AUDIT"
    elif (stalled or slow_momentum) and max(blocks)==blocks[0] and Path("bin/momentum/finish.epoch").is_file():
        decision="COMPARE_SCALAR_JOINT_WITH_MOMENTUM_AT_SAME_BASE"
        report["trigger"]="terminal scalar direction limit" if stalled else "slow momentum after full scalar window"
    elif terminal and not (stalled and max(blocks)==blocks[0]):
        decision="SCALAR_TERMINATION_DID_NOT_IDENTIFY_MOMENTUM_LIMIT"
elif terminal:
    decision="NO_ACCEPTED_FIELD"
if decision is not None:
    report["decision"]=decision
    save(Path("momentum-decision.json"),report)
R11_MOMENTUM_PICK
    if ! test -s momentum-decision.json; then sleep 30; fi
done
seed="$(python3 - <<'R11_MOMENTUM_SEED'
import json
from pathlib import Path

record=json.loads(Path("momentum-decision.json").read_text())
if record["decision"]=="COMPARE_SCALAR_JOINT_WITH_MOMENTUM_AT_SAME_BASE":
    print(record["packet"]["path"])
R11_MOMENTUM_SEED
)"
if test -n "$seed"; then
    remaining=$((10800-$(date +%s)+$(cat start.epoch)))
    cooperative=$((remaining-1850))
    status=0
    bash solve.sh momentum momentum "$seed" momentum Reanchor \
        "$cooperative" 12 .00625 "$r/frame/recycled.json" 12 || status=$?
    case "$status" in
        0|1|2) test -s momentum/state.steady.json ;;
        *) exit "$status" ;;
    esac
fi
date +%s > momentum.finish.epoch
