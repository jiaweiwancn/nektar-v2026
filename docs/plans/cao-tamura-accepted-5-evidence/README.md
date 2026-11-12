# CT-Accepted-5 evidence

Start with `contract.json`, `qualified-v1/final-decision.json` and the
[package report](../cao-tamura-accepted-5.md). `CT4-observation-correction.md`
supplements CT4 without editing its sealed records.

The Linux experiment root is
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-5-20261112` (R5).
CT4 data/runtime remain under the sibling `ct-accepted-4-20261112` (R4).
Native source is pinned to `123f31c974e0b75e551e43a5203749a89378fade`;
final CT5 collection source is `acafc765f24392446a2fea0c4a66fd458152d1bb`.
The separately tested package-wide half-step guard is committed at
`f8473f0`; numerical classification functions are AST-identical. Actual
per-request source and runtime hashes take precedence over these labels.

## Contents

- `qualified-v1/compact-v1/*-complete-comparisons.json`: complete model/case
  initialization and milestone bodies, original ambiguous A93 pair, actual
  signed traction profiles, timestep gates, fields, domain witnesses, budget
  primitives, consumed boundaries and state/native-history identities.
- `qualified-v1/compact-v1/load-and-stage-diagnosis.json`: all five original
  accepted states, term reconstruction including open loads, signed radial
  localization, analytical weak balances and actual six-stage traction
  decomposition. Full coefficients/basis packets remain at R5/diagnosis-v1.
- `qualified-v1/post-v1/results.json`: common-grid equal-time initial/evolved
  contrasts, stage scalar inventories, independent rough pressure actions,
  and reconstructed outlet mismatch distinguished from natural data.
- `qualified-v1/compact-v1/observation-correction.json`: old/new 40 nominal
  distances, 80 envelope changes, all 36 Fig.7 ratios, smooth classification
  flips, independent geometry controls and separate inner/outer scale ratios.
  The full 796-point audit is sealed at R5/observations-v1/results.json; its
  hash is included.
- `qualified-v1/compact-v1/negative-controls.json`: intended failures and
  discriminating reasons. Default replay byte-equivalence and side controls
  are also in `qualified-v1/execution-and-preparation/preflight-controls-v1.json`.
- `qualified-v1/compact-v1/provenance.json`, `actual-loaded-libraries.json`:
  source/spec/data identities, exact immutable requests, binaries, actual
  libraries, build flags and complete heavy-artifact bindings.
- `qualified-v1/compact-v1/native-costs.json`: exclusive instrumented scopes
  and native process CPU/wall/RSS for 13 distinct evolution segments.
- `qualified-v1/completed-resource-receipts.json`: closed numerical/diagnostic
  launch receipts. `resource-highwater-v1.json` adds time/PID/thread ownership
  and suspension records. `archive-resource-final.json` and
  `highwater-export-resource-final.json` close the two final bookkeeping
  families. The collector's pending-own-finalizer snapshot is historical;
  these completed receipts close it.
- `qualified-v1/execution-and-preparation`: initial audit and exact recorded
  command bodies, source synchronization, preserved preparation failures and
  recovery exit statuses. Process completion is distinct from qualification.
- `halfstep-limit-v1.json` and its resource/command receipts: the global
  conditional-control cap, distinct from the trajectory-producing versions.
- `linux-export-sha256.json`, `transfer-verification.json`: original Linux
  compact archive identity and verified nested request/artifact seals. All
  original exported JSON bytes are retained. `SHA256SUMS.json` additionally
  covers the whole committed evidence tree except itself.

No matrices, coefficient vectors, checkpoints or dense basis arrays are
copied into Git. Hashes and exact Linux paths locate them. The initial member
receipts preserve native counts, per-block counts, configuration and reference
hashes, raw/effective state receipts and both live boundary hashes. Same-space
timestep initial coefficients agree; cross-space physical functions agree
while their projections may differ.

## Reproduction

Use Linux and a new controlled destination for every new job. Do not rerun
inside a sealed output directory. The exact commands actually executed are in
the launch receipt `execution.command` arrays, including their Python command
bodies. These bind the historical source versions used for explicit recovery.
The following read-only rechecks use the retained complete raw data and final
source. Replace `UNIQUE` before running; all numerical libraries must be one
thread **before** launching Python.

```bash
R=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-5-20261112
C=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-4-20261112
T="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
RT="$C/build-v3/runtime"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$R/tmp" LD_LIBRARY_PATH="$RT/lib"

python3 -B "$T/ct_accepted3_command.py" --control "$R/launch-post-UNIQUE" -- \
  python3 -B "$T/ct_accepted5_post.py" --root "$R" --ct4 "$C" \
  --output "$R/post-UNIQUE"
python3 -B "$T/ct_accepted3_command.py" --control "$R/launch-collect-UNIQUE" -- \
  python3 -B "$T/collect_ct_accepted_5.py" --root "$R" --ct4 "$C" \
  --runtime "$RT" --campaign "$R/controller-halfstep-v1" \
  --output "$R/compact-UNIQUE"
```

For new evolution, create a fresh isolated experiment/source root, synchronize
the intended commit, and run `ct_accepted5_observations.py`,
`ct_accepted5_diagnose.py`, then `ct_accepted5_run.py` with that root, unchanged
R4 references/runtime and `contract.json`. Each CLI documents its required
paths. Keep the one-split selection and single conditional .0005 rule.
`ct_accepted5_halfstep.py` is the explicit, tested conditional A93 control;
its request binds the original ambiguous pair and refuses a second quarter
initialization. The executed recovery used `ct_accepted5_resume.py` only to
reuse verified completed members after reporting errors. It did not restart
an old physical failure. Native `.02` segments consume their own unchanged-dt
`.008` histories; timestep changes always have independent zero-step histories.

Kernel/specification/data reuse is verified against CT4 source hashes, and
the native C++ diff from the reviewed baseline is empty. CT1–CT4 numerical
and resource failures, failed preparations and the historical pressure flag
are not relabelled. Native 6 mm bounded accuracy is distinct from physical
surface adequacy, smooth preparation, precursor admission and hill validation.
