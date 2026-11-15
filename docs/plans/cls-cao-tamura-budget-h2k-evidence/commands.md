# H2k execution and replay commands

16 November 2026, Asia/Shanghai. Numerical work used Linux. These are immutable
completed paths; use new unused output directories for any reproduction. The
machine command/packet index is `commands.json`; `artifact-index.json` binds every
retained numerical input/output. No long continuation is part of these commands.

```bash
R='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2k-20261115'
J='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2j-20261115'
G='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2g-20261113'
H1='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111/compatibility-resolved'
T="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
PROBE="$R/probe-runtime/RANSWallResolutionProbe"
SOLVER="$G/runtime-v2/IncNavierStokesSolver"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_THREAD_LIMIT=1
export MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export TMPDIR="$R/tmp"
findmnt -no SOURCE,UUID,TARGET --target "$R"
# Required UUID: 5AB8A669B8A642FD. Account for other project workers first.
```

Source came from a symlink-preserving isolated copy of H2j, with the exact Mac
changed blobs overlaid and checked. A new build avoids copied-mtime stale objects.
The actual configuration/build was:

```bash
cmake -S "$T" -B "$R/build" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS' \
  '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake' \
  -DRANS_BUILD_SOLVER=OFF \
  '-DRANS_PREBUILT_SPATIAL_LIBRARY=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2b-20261112/curved-lazy-final/runtime/lib/libSpatialDomains.so.5.9.0'
cmake --build "$R/build" --parallel 4 --target RANSWallResolutionProbe
```

`continue_cls_epawf_pimple.freeze()` copied the new probe and linked closure to
`probe-runtime` before qualification. The deployed physical solver was reused
unchanged. Full identities are in `source-runtime.json`, the certificate and the
indexed `probe-runtime/identity.json`. Do not rebuild or refreeze over these paths.

The adaptive driver's phases are `patch`, `construct`, `context`, `search` and
`reserve`. Shared arguments are `--root "$R" --parent "$G" --previous "$J"`;
native phases also use `--probe "$PROBE"`. The retained records specify:

| Operation | Inputs / controls | Immutable output |
|---|---|---|
| Upper patch | geometries60–77, H2j actual dual witnesses | patch1 |
| Exact strong-null attempt | no independent retained complement | construct1 |
| Compensated construction | `--patch "$R/patch1" --compensated --weights .01,1` | construct2 |
| Signed native context check | `--construction "$R/construct2"` | epoch1-context |
| Anchor A then B, buffered/actual, independent replay | 78-column context | epoch1-search |
| Outer patch | geometries36–77, `--retained-context "$R/epoch1-context/context.json" --proof-source "$R/epoch1-search" --target-count 16` | patch2 |
| Targeted construction | `--patch "$R/patch2" --retained-context "$R/epoch1-context/context.json" --compensated --weights 1 --divergence-compensator` | construct3 |
| Signed native context check | construct3 and retained epoch1-context | epoch2-context |
| Unscaled LP, equivalent scaled LP, exact-cap/SOC organization | same frozen95-column context; all attempts charged | epoch2-search* |
| Successful complete separation | `--exact-caps --equilibrated`, saved initial coefficients | epoch2-search-SOC |
| One stricter interior construction | `reserve --context "$R/epoch2-context/context.json" --candidate "$R/epoch2-search-SOC/candidate.json" --exact-caps` | startup-reserve |

All commands include `--output <unused output>`. `result.json`, event files,
terminal coefficients and executed-source hashes preserve failed and successful
numerical organizations. The interrupted scaled LP is explicitly charged and
recorded in `epoch2-stop-for-exact-quadratics.json`; it is not an obstruction.
The final exact physical seed is independently reproducible from `selected/seed.json`
and its hash, without depending on optimizer bitwise reproducibility across SciPy
versions. A new search uses its own pre-build wall origin and aggregate ledger;
the certificate-bound completed ledger must not be edited or reused as new work.

After source/runtime freeze, the actual qualification and startup sequence was:

```bash
cd "$T"
python3 qualify_cao_tamura_budget_h2k.py seen --root "$R" \
  --context "$R/epoch2-context/context.json" --candidate "$R/startup-reserve/candidate.json" \
  --probe "$PROBE" --preflight "$G/comparison.json" --h1 "$H1" \
  --output "$R/seen-reserve"
python3 qualify_cao_tamura_budget_h2k.py fresh-import --root "$R" \
  --context "$R/epoch2-context/context.json" --candidate "$R/startup-reserve/candidate.json" \
  --seen "$R/seen-reserve/receipt.json" --probe "$PROBE" --solver "$SOLVER" \
  --preflight "$G/comparison.json" --h1 "$H1" \
  --source-revision c47630090bcb4d52a86c8ba5eeb5976a9ade7264 \
  --output "$R/fresh-epoch2"
# Checkpoint2 was committed/pushed before this physical invocation.
python3 qualify_cao_tamura_budget_h2k.py pilot --root "$R" \
  --context "$R/epoch2-context/context.json" \
  --certificate "$R/fresh-epoch2/certificate.json" --solver "$SOLVER" \
  --preflight "$G/comparison.json" --h1 "$H1" --output "$R/pilot"
```

The last command's ordinary coarse/fine/split trajectories all pass. Its original
wrapper stops during the expected deliberate rejection because of the diagnosed
empty-local-checkpoint validation defect. That ledger remains intact. The
validation-only repair is tested and requalified with **no ordinary trajectory
rerun and no fresh-grid consumption**:

```bash
python3 test_cao_tamura_budget_h2i.py --output "$R/rollback-validator-regression.json"
python3 revalidate_cao_tamura_budget_h2k.py --root "$R" --solver "$SOLVER" \
  --certificate "$R/fresh-epoch2/certificate.json" --preflight "$G/comparison.json" \
  --h1 "$H1" --output "$R/pilot-requalified"
python3 analyse_cao_tamura_budget_h2k.py --root "$R" --probe "$PROBE" \
  --comparison "$G/comparison.json" --output "$R/post-startup-final"
python3 collect_cao_tamura_budget_h2k.py --root "$R" \
  --comparison "$G/comparison.json" --output "$R/final-evidence-v2"
```

The repaired validator explicitly uses the accepted parent for a deliberate
rejection, while ordinary invocations still require their actual accepted local
checkpoint. It additionally checks signature, complete histories/time labels,
stage data, diagnostic/sample bits and unchanged parent bytes. Its new focused
rollback is the second deliberate rejected attempt. The two failed analysis
setup attempts concern only endpoint last-bit time labels and absent zero-step
reporter samples; neither ran a native solve. The completed strict readback uses
native h0 samples and preserves actual time labels.

Focused native-context and mathematical tests are indexed in `verification.json`.
`test_cao_tamura_budget_h2k.py` takes `--context`, `--previous`, `--output`;
the LP, exact-cap, handoff-routing and affected H2i tests take `--output`.
GL28 was consumed exactly once for the frozen hill seed. Replaying its archived
packet is a repeatability check, never another fresh qualification.

GNU-time records use `%e %U %S %M` (RSS in KiB) or `-v`. Native JSONL timers retain
inclusive and exclusive costs separately. Detailed failed/selected attempt packets,
physical checkpoints, native histories and rejection evidence remain in Linux
storage. No write/delete policy or production hard audit was weakened.
