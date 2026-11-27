#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
old="$(dirname "$r")/wj-c1-b3-bounded.xjpRylk7"
export PYTHONPATH="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
while ! test -s phase-one.finish.epoch; do
    test "$(($(date +%s)-$(cat start.epoch)))" -lt 5200
    sleep 15
done
while ! bash capacity.sh phase-two-selection.ps 2; do
    test "$(($(date +%s)-$(cat start.epoch)))" -lt 5200
    sleep 15
done
seed="$(python3 - "$old/bounded-active/state.accepted-000003.json" <<'R11_PICK'
from pathlib import Path
import sys
from c1_bounded_trials import load, packed, save
from resources import identity

start=Path(sys.argv[1])
choices=[start,*Path("retained").glob("state.accepted-*.json")]
values=[max(load(f)["original_residuals"]) for f in choices]
base=values[0]
best, seed=min(zip(values,choices))
assert best>1.e-10, "retained span reached native root tolerance; qualify it first"
assert (base-best)/base<.01, "retained span shows useful progress; reassess before assembling an inverse"
packet=load(seed)
coef=packet["state"]["state_sha256"]
packed(packet,coef)
selection=Path("matched-seed.json")
record={"packet":identity(seed.resolve()),"coefficients":coef,"native_max":best,
        "relative_decrease_in_retained_span":(base-best)/base}
if selection.exists():
    assert load(selection)==record
else:
    save(selection,record)
print(seed.resolve())
R11_PICK
)"
if test -s retained/state.accepted-000001.json && ! test -s resume-controls.json; then
    bash resume-controls.sh >> resume-controls.log 2>&1
fi
if ! test -s bin/coupled/finish.epoch; then
    bash build.sh coupled 2 > build-coupled.log 2>&1
fi
remaining=$((10800-$(date +%s)+$(cat start.epoch)))
cooperative=$((remaining-2190))
test "$cooperative" -gt 2300
status=0
bash solve.sh coupled scalars "$seed" scalar-blocks Reanchor \
    "$cooperative" 21 .00625 "$r/frame/recycled.json" 12 || status=$?
test "$status" = 0 || { test "$status" = 2 && test -s scalars/state.steady.json; }
date +%s > phase-two.finish.epoch
