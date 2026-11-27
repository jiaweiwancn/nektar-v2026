#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
old="$(dirname "$r")/wj-c1-b3-bounded.xjpRylk7"
while ! test -s bin/ratio/finish.epoch; do
    test "$(($(date +%s)-$(cat start.epoch)))" -lt 3500
    sleep 15
done
test -s reduced/finish.epoch
bash replay.sh ratio terminal --frame-output "$r/frame"
bash replay.sh ratio refined --proposal "$r/reduced/legacy-025" \
    --proposal "$r/reduced/ratio-.025" --proposal "$r/reduced/ratio-.00625" \
    --proposal "$r/reduced/ratio-.003125"
status=0
bash solve.sh ratio retained "$old/bounded-active/state.accepted-000003.json" retained \
    Reanchor 1500 6 .00625 "$r/frame/recycled.json" 12 || status=$?
test "$status" = 0 || { test "$status" = 2 && test -s retained/state.steady.json; }
date +%s > phase-one.finish.epoch
