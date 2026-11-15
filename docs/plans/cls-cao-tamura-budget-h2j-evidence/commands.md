# H2j executed commands

15 November 2026, Asia/Shanghai. These paths identify completed immutable outputs.
For reproduction, choose a new unused result root and synchronize the exact
implementation blobs in `commit-binding.json`; never overwrite these results.
Numerical conclusions and timings come from Linux. All runs use one rank and
one computational-library thread. No production solver was rebuilt.

```bash
R='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2j-20261115'
I='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2i-20261115'
G='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2g-20261113'
H1='/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111/compatibility-resolved'
T="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
PROBE="$R/probe-runtime/RANSWallResolutionProbe"
SOLVER="$G/runtime-v2/IncNavierStokesSolver"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export TMPDIR="$R/tmp"
findmnt -no SOURCE,TARGET,UUID /media/jay/新加卷1
# Required UUID: 5AB8A669B8A642FD. Inspect other workers before scheduling.
```

The H2i source tree was copied with symlinks preserved, then overlaid with the
Mac changes after checksum verification. A fresh build prevents stale objects
from copied mtimes. Configuration and the successful build:

```bash
cmake -S "$T" -B "$R/build" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS' \
  '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake' \
  -DRANS_BUILD_SOLVER=OFF \
  '-DRANS_PREBUILT_SPATIAL_LIBRARY=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2b-20261112/curved-lazy-final/runtime/lib/libSpatialDomains.so.5.9.0'
cmake --build "$R/build" --parallel 4 --target RANSWallResolutionProbe
```

The first build caught misleading indentation under `-Werror`; the corrected
build passed. `continue_cls_epawf_pimple.freeze()` copied the probe and linked
closure into `probe-runtime` before numerical qualification. Its original and
frozen hashes are in `probe-runtime/identity.json`; selected hashes are in
`source-revisions.json`. Historical solver and library files were preserved.

Executed sequence, with independent preparation, flat and handoff jobs allowed
to run concurrently under the stated resource controls:

```bash
cd "$T"
python3 run_cao_tamura_budget_h2j.py prepare --parent "$G" --previous "$I" \
  --probe "$PROBE" --output "$R/prepare"
python3 run_cao_tamura_budget_h2j.py handoff --parent "$G" \
  --probe "$PROBE" --solver "$SOLVER" --output "$R/handoff-final"
python3 run_cao_tamura_budget_h2j.py flat --parent "$G" \
  --probe "$PROBE" --output "$R/flat"
python3 run_cao_tamura_budget_h2i.py search --parent "$G" \
  --context "$R/prepare/level1/context.json" --output "$R/joint-search"
for anchor in A B; do
  python3 run_cao_tamura_budget_h2i.py repeatable \
    --candidate "$R/joint-search/$anchor-actual/best-candidate.json" \
    --context "$R/prepare/level1/context.json" --probe "$PROBE" \
    --preflight "$G/comparison.json" --h1 "$H1" --output "$R/native-$anchor"
done
python3 analyse_cao_tamura_budget_h2j.py --root "$R" --previous "$I" \
  --output "$R/compact"
python3 test_cao_tamura_budget_h2j.py \
  --receipt "$R/handoff-final/issuance-receipt.json" \
  --output "$R/positive-handoff-qualified-tests.json"
python3 test_cao_tamura_budget_h2i.py --output "$R/h2i-regression.json"
```

The search uses the original H2i scalar retention, fixed-k anchors and eight
epsilon directions, plus the frozen 44-column velocity basis. Actual-limit
problems remove buffered cuts/bounds. The reducer independently reconstructs
necessary cuts and conservative dual lower bounds in the augmented space; six
replay/diagnostic scans are charged to the aggregate. This is not another
optimizer or a hill attempt. `optimizer.json` records search-internal native
checks as zero; the two separate native calls are recorded in `native.json`.

Timing wrappers use `/usr/bin/time -f '%e %U %S %M' -o <usage.txt>`. They include
driver setup, native processing and output. The process RSS unit is KiB. The
initial `handoff` native import succeeded but its Python wrapper rejected an
empty success string. The corrected driver and producer were rerun cleanly at
`handoff-final`; both outputs and the first log remain indexed.

After reduction, the independent witness helper was corrected to retain full
`-R:G`, including isotropic divergence work. Compact witness contractions were
updated from the retained tensors without another optimizer or inventory scan.
Native governing production, budgets and tensor data were already correct.

Final collection verifies exact source hashes against the Mac implementation
commit, preserves full runtime closures, and hashes the native artifact tree.
The lightweight `strict-zero-enrichment-replay.json` independently compares all
old response blocks, memberships and weights as binary64 bytes. All selected
source bytes, including signed-zero handling, match the committed implementation.

## Conditional operations not executed

No candidate passes complete covariance qualification. Therefore `fresh` was not
called: the single hill GL28 audit, hill certificate/import, coarse8/fine16,
coarse4+4 and curved rollback remain unexecuted. The separately qualified flat
V2 fixture has its own base and ledger. It does not reserve or consume hill GL28.
No long continuation was launched. A failed phase-I field must not be passed to
the startup runner or treated as an accepted physical checkpoint.
