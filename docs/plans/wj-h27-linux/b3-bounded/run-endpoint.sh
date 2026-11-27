#!/usr/bin/env bash
set -euo pipefail
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
stage=$1
seed="$(realpath "$2")"
name=$3
mode=$4
seconds=$5
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite:$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h26"
table="/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/source-v2/docs/plans/wj-h27-linux/task2r-r1/tables-exact-inner-v1.json"
paired="/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-paired-axial.UXeQ6Noe"
q="$r/$name"
bash workers.sh "$name-workers.ps" 66
case "$mode" in
  pending|pending-all|qualify)
    printf '%s  %s\n' 3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1 "$table" | sha256sum -c -
    /usr/bin/time -v -o "$name.setup.cost" python3 "$suite/c1_steady_case.py" \
      "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/hybrid-fine/run.xml" "$q" --seed "$seed" --eta 0 --quadrature 16 24 \
      --seconds 190 --evaluate-only --check-actions
    bash "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" \
      "$q" "$r/bin/$stage/IncNavierStokesSolver" 280
    packet="$q/state.diagnostic-steady.json"
    python3 - "$seed" "$packet" "$mode" <<'CHECK_RICH_R10'
from c1_bounded_trials import SETTINGS,load,packed
from pathlib import Path
import sys
a,b=(load(Path(p)) for p in sys.argv[1:3])
packed(a,a["state"]["state_sha256"])
packed(b,a["state"]["state_sha256"])
assert all(a["state"][key]==b["state"][key] for key in SETTINGS)
assert (b["state"]["quadratureX"],b["state"]["quadratureR"])==(16,24)
print("coefficient",b["state"]["state_sha256"],"original norm",max(b["original_residuals"]))
if sys.argv[3]=="qualify":
    assert max(a["original_residuals"])<=1e-10 and max(b["original_residuals"])<=1e-10
CHECK_RICH_R10
    if [ "$mode" = pending ]; then
      /usr/bin/time -v -o "$q/equations.cost" timeout --kill-after=15s "$((seconds+25))" \
        python3 "$suite/c1_oracle.py" --diagnostic --native "$packet" --table "$table" --order 28 \
          --out "$q/equations.json" --seconds "$seconds"
      exit 0
    fi
    /usr/bin/time -v -o "$q/equations.cost" timeout --kill-after=15s "$((seconds+25))" \
      python3 "$suite/c1_oracle.py" --diagnostic --native "$packet" --table "$table" --order 28 \
        --out "$q/equations.json" --seconds "$seconds" &
    eq=$!
    /usr/bin/time -v -o "$q/margins.cost" timeout --kill-after=15s "$((seconds+25))" \
      python3 "$suite/c1_margins.py" --diagnostic --native "$packet" --table "$table" --order 28 \
        --one-sided-faces --out "$q/margins.json" --seconds "$seconds" &
    mg=$!
    eq_status=0
    mg_status=0
    wait "$eq" || eq_status=$?
    wait "$mg" || mg_status=$?
    test "$eq_status" = 0 && test "$mg_status" = 0
    /usr/bin/time -v -o "$q/points.cost" timeout --kill-after=15s 310 \
      python3 "$suite/c1_refined_witness.py" inspect --native "$packet" --table "$table" --order 28 \
        --seconds 280 --out "$q/points.json" --release-controls --minima inner "$q/margins.json" \
        --minima outer "$q/margins.json" --minima faces_inner "$q/margins.json" \
        --minima faces_outer "$q/margins.json" --minima forcing_support "$q/margins.json"
    python3 - "$q" "$table" <<'B3_R10_ADMISSION' > "$q/admission.txt"
from pathlib import Path
import sys
from c1_steady_continue import original_admission
q,table=map(Path,sys.argv[1:])
print(original_admission(q/"state.diagnostic-steady.json",table,
                         {n:q/(n+".json") for n in ("equations","margins","points")}))
B3_R10_ADMISSION
    cat "$q/admission.txt"
    ;;
  paired)
    test "$(cat "$q/admission.txt")" = QUALIFIED
    printf '%s  %s\n' 73c264cf3d237311049fec71f07c32fa44a5b38c2cd2cd308c85103c9f62ab98 \
      "$paired/frozen-core.json" | sha256sum -c -
    args=()
    for f in A2 B2 A3; do
      p="$paired/$f-rich"
      args+=(--field "$f" "$p/state.diagnostic-steady.json" --margins "$f" "$p/margins.json")
    done
    args+=(--field B3 "$q/state.diagnostic-steady.json" --margins B3 "$q/margins.json")
    /usr/bin/time -v -o paired-level3.cost timeout --kill-after=15s "$((seconds+25))" \
      python3 "$suite/c1_paired_profiles.py" --frozen "$paired/frozen-core.json" --table "$table" \
        --out "$r/paired-level3" --seconds "$seconds" --coverage --noise-receipt "$paired/baseline-profile/paired-comparison.json" \
        --retained-patches "$paired/level-2-profile/paired-comparison.json" \
        --retained-patches "$paired/partial-third-profile/paired-comparison.json" \
        --patch -6.25 .15 .09 .012 --patch .17578125 .10041061001735275 .04 .008 \
        --compare A2 A3 --compare B2 B3 --compare A3 B3 "${args[@]}"
    ;;
  *) exit 2 ;;
esac
