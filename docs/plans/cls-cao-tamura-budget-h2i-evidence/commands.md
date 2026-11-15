# H2i reproduction commands

Recorded 15 November 2026, Asia/Shanghai. Use an unused output root for a replay;
the paths below describe the completed artifacts and must not be overwritten.
The tested source commit and individual source hashes are in `commit-binding.json`
and `source-revisions.json`. Numerical runs are Linux-only. SSH used the existing
Tailscale fallback when `linux-usb` was unavailable.

```bash
R='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2i-20261115'
G='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2g-20261113'
H1='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111/compatibility-resolved'
T="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
PROBE="$R/probe-runtime/RANSWallResolutionProbe"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export TMPDIR="$R/tmp"
findmnt -no SOURCE,TARGET,UUID /media/jay/新加卷1
# Required UUID: 5AB8A669B8A642FD. Check active workers before a new build.
```

Source was overlaid from the Mac onto a separate copied H2h source tree with
matching checksums, then configured in a **fresh** build directory. This avoids
stale objects from copied timestamps. The external test project uses installed
Nektar++ dependencies and the already verified SpatialDomains closure:

```bash
cmake -S "$T" -B "$R/build" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS' \
  '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake' \
  -DRANS_BUILD_SOLVER=OFF \
  '-DRANS_PREBUILT_SPATIAL_LIBRARY=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2b-20261112/curved-lazy-final/runtime/lib/libSpatialDomains.so.5.9.0'
cmake --build "$R/build" --parallel 4 --target RANSWallResolutionProbe
```

`continue_cls_epawf_pimple.freeze()` copied the probe and linked closure to
`probe-runtime`; both original and frozen identities are retained. Do not
refreeze or rebuild an executable underneath a certificate or running process.
The unchanged production solver is `$G/runtime-v2/IncNavierStokesSolver`.

Executed preparation/search/check sequence (names are original durable outputs):

```bash
cd "$T"
python3 run_cao_tamura_budget_h2i.py prepare --parent "$G" --probe "$PROBE" \
  --h1 "$H1" --preflight "$G/comparison.json" --output "$R/joint-context"
python3 run_cao_tamura_budget_h2i.py search --parent "$G" \
  --context "$R/joint-context/context.json" --output "$R/joint-search"
for anchor in A B; do
  python3 run_cao_tamura_budget_h2i.py repeatable \
    --candidate "$R/joint-search/$anchor-actual/best-candidate.json" \
    --context "$R/joint-context/context.json" --probe "$PROBE" \
    --preflight "$G/comparison.json" --h1 "$H1" --output "$R/native-$anchor"
done
python3 run_cao_tamura_budget_h2i.py flat --probe "$R/build/RANSWallResolutionProbe" \
  --preflight "$G/comparison.json" --h1 "$H1" --output "$R/flat-H1"
python3 analyse_cao_tamura_budget_h2i.py --root "$R"
```

The flat readback's build-path executable and dependency hashes equal the frozen
probe's content hashes. Its actual original paths are recorded. Preparation,
search, candidate checks and flat readbacks were timed with
`/usr/bin/time -f '%e %U %S %M' -o <usage.txt>`; native receipts also retain wall
time. The independent reducer was run twice, with cumulative scan accounting,
to add explicit matching-side diagnostic selection. It solves no LP.

The reducer saves two bounded accessibility packets. Each was passed to
`diagnose_cls_accepted_covariance.run_probe()` using the corresponding native
check's exact XML and its SHA256, with output roots `accessibility-A/B`.
`collect_cao_tamura_budget_h2i.py` then contracts returned responses using the
saved physical mass/divergence forms. This does not enrich the basis or alter
the candidate. Original input/output packets, checksums and commands are in the
artifact index.

Focused checks executed on Linux:

```bash
python3 test_cao_tamura_budget_h2i.py --output "$R/unit-qualified.json"
python3 test_cao_tamura_h2h_scalar.py --output "$R/h2h-math-regression.json"
cmake --build "$R/build" --parallel 4 --target RANSCaoTamuraBudgetH2gChecks
"$R/build/RANSCaoTamuraBudgetH2gChecks" \
  '/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2c-20261112/preflight.json' \
  "$R/calibration-regression.json"
python3 collect_cao_tamura_budget_h2i.py --root "$R" --parent "$G" \
  --output "$R/compact-evidence" --revision '<tested commit or exact source hash manifest>'
```

Final collection verifies source/runtime hashes and measures typed-array hashing
I/O separately. Full physical samples and raw candidates stay on Linux; compact
evidence references their checksummed native packets. The complete scans used by
independent reductions are charged separately from optimizer scans.

## Conditional operations, deliberately not executed

`fresh` requires a frozen, seen/native-qualified candidate, matched context,
basis and runtime, then reserves one GL28 ledger entry before invoking the probe.
`pilot` requires the passing exact GL28 receipt and a model-specific V2
certificate. Failed or mismatched receipts cannot launch it. Neither prerequisite
was met in H2i. Do not run GL28 to tune another candidate, edit a checkpoint's dt,
or promote either phase-I diagnostic state to an accepted seed.

The qualified path checks coarse h0→coarse8→fine h0→fine16→coarse4→restart4→
deliberate rollback **after each invocation**. Any unexpected status, report,
rejection, equation count, solved-epsilon or checkpoint mismatch is terminal.
No long continuation is authorized by this package.
