#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
export PYTHONPATH="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
bash workers.sh model-refinement-workers.ps 67
python3 refine_model_candidate.py bounded-enriched "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/hybrid-fine/run.xml" model-refinement
bash "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" "$r/model-refinement" "$r/bin/recycled/IncNavierStokesSolver" 425
