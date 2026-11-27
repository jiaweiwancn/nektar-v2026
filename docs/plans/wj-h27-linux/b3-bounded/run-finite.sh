#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
suite="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
export PYTHONPATH="$suite"
test -s bin/finite/finish.epoch
bash workers.sh finite-workers.ps 68
sha256sum -c archive.sha256
python3 c1_bounded_trials.py prepare "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM" "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-frozen-size.pnvR6Uif" "$r/frozen-trials" --seconds 750
date +%s > frozen-trials/start.epoch
bash "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" "$r/frozen-trials" "$r/bin/finite/IncNavierStokesSolver" 850
python3 c1_bounded_trials.py summarize "$r/frozen-trials"
date +%s > frozen-trials/finish.epoch
