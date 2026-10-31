#!/usr/bin/env bash
# Read-only progress for a CBC100 preflight or study, from the Mac or Linux.
set -euo pipefail
hit_run=${1:?Usage: check_band100_status.sh /media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/RUN}
hit_tool=/media/jay/新加卷/nektar-v2026/solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/run_band100.py
if [[ $(uname -s) == Linux && -f "$hit_tool" ]]; then
    exec python3 "$hit_tool" --status "$hit_run"
fi
printf -v hit_quoted_run '%q' "$hit_run"
exec ssh -o BatchMode=yes -o ConnectTimeout=15 jay@100.105.247.110 \
    "python3 '$hit_tool' --status $hit_quoted_run"
