#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
stage=$1
if test -s momentum-intent.txt; then
    while ! test -s momentum.finish.epoch; do
        test "$(($(date +%s)-$(cat start.epoch)))" -lt 10010
        sleep 15
    done
fi
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
old="$(dirname "$r")/wj-c1-b3-bounded.xjpRylk7"
readarray -t endpoint < <(python3 - "$old/bounded-active/state.accepted-000003.json" <<'R11_ENDPOINT'
import sys
from pathlib import Path
from c1_bounded_trials import load, packed, save
from resources import identity

old=Path(sys.argv[1])
choices=[old,*Path("retained").glob("state.accepted-*.json"),*Path("scalars").glob("state.accepted-*.json"),
         *Path("momentum").glob("state.accepted-*.json")]
best=min(choices,key=lambda p:max(load(p)["original_residuals"]))
packet=load(best)
norm=max(packet["original_residuals"])
coefficients=packet["state"]["state_sha256"]
packed(packet,coefficients)
mode="qualify" if norm<=1e-10 else "reused-original-base" if best==old else "pending"
save(Path("endpoint-selection.json"),{"packet":identity(best.resolve()),"coefficient_sha256":coefficients,
                                     "native_max":norm,"mode":mode})
print(best.resolve())
print(mode)
R11_ENDPOINT
)
test "${#endpoint[@]}" = 2
mode="${endpoint[1]}"
if [[ "${endpoint[0]}" == "$r/momentum/"* ]]; then stage=momentum; fi
if test "$mode" = reused-original-base; then
    exit 0
fi
driver="$r/source/docs/plans/wj-h27-linux/b3-bounded/run-endpoint.sh"
sed 's/bash workers.sh "$name-workers.ps" 66/bash capacity.sh "$name-workers.ps" 3/' "$driver" > endpoint-driver.sh
remaining=$((10800-$(date +%s)+$(cat start.epoch)))
if test "$mode" = qualify; then
    test "$remaining" -gt 1630
    seconds=1370
else
    test "$remaining" -gt 715
    seconds=$((remaining-230))
    if test "$seconds" -gt 1000; then seconds=1000; fi
fi
bash endpoint-driver.sh "$stage" "${endpoint[0]}" endpoint-rich "$mode" "$seconds"
date +%s > endpoint-rich/finish.epoch
remaining=$((10800-$(date +%s)+$(cat start.epoch)))
if test "$mode" = qualify && test "$(cat endpoint-rich/admission.txt)" = QUALIFIED && test "$remaining" -gt 800; then
    bash endpoint-driver.sh "$stage" "${endpoint[0]}" endpoint-rich paired 620
fi
