# H18 checkpoint-audit reproduction

The authoritative completed root is
`/home/jay/build/wj-h18-checkpoint-audits-20261114`. Its ledger is CLOSED.
These commands describe the executed package; they do not authorize reopening
that root or another scientific trajectory. Keep heavy inputs in place.
Source/Git work uses the Mac WJ worktree; Linux builds/numerics use the isolated
capsule, GNU 11.4 Release and the qualified private loader/library closure.

## Configuration and roles

Clone the manifest-selected H14 session; change execution controls only:

```xml
<I PROPERTY="RANSWallRetention" VALUE="BoundedV1" />
<I PROPERTY="RANSWallMonitoringPolicy" VALUE="CheckpointAuditV1" />
<I PROPERTY="RANSWallEvaluationReuse" VALUE="ExactStepV1" />
<I PROPERTY="RANSWallImmutableRowReuse" VALUE="Off" />
<I PROPERTY="RANSWallPimpleSnapshots" VALUE="Off" />
<I PROPERTY="RANSWallMonitorIntent" VALUE="/owned/run/intent.json" />
```

Diagnostics and StopOnBound stay On, PimpleVerify stays Off. The immutable
`WJMonitorIntentV1` includes the existing nonce/source/observer/profile/qh0 pins,
plus:

```json
{
  "recovery_stride": 1,
  "scientific_stride": 4,
  "maximum_absolute_step": 4,
  "dt_bits": "3f647ae147ae147b",
  "substep_evidence_steps": [],
  "preserve_qualification_artifacts": true
}
```

The shown dt bits are .0025. Strides use absolute steps, not local invocation
indices. Scientific checkpoints are also recovery checkpoints. Final state and
restored input are always audited. Optional `terminal_start`/`terminal_stride`
refine the scientific cadence. Recovery10/scientific100/terminal50 are future
sampling proposals, not acceptance criteria or duration admission.

`substep_evidence_steps` selects rich k evidence; bound events capture the
operator while available. It never disables the bound. Late discarded operators
are not recreated. The H18 pilot preserves every recovery; `rotate_owned`
accepts only pinned owned recoveries explicitly superseded and not protected,
retaining newest 3. Pending, scientific, verified and first-failure anchors remain.
Deletion/rotation is tested with synthetic files, not thousands of CFD steps.

A recovery has State.bin, WJCompletedCallbackV1 and WJNativeCompletionV1 pins.
The latter binds actual state/history/time/alpha, reference/defect/original qh0,
producer/configuration and runtime/source identity. Cursors distinguish native
completion, durable recovery and committed verification. Pending strict requests
use `import_role=PendingNativeV1` plus `native_completion`; qualified requests
continue to require their native verification commit. Restored completed states
use their producing callback, not a zero-alpha CP0 audit. Startup verification
has its own commit and never rewrites original audit status.

## Recorded commands

All six numerical-library thread variables were 1. Compile parallelism was 8;
host total stayed below 72. The focused build used:

```sh
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
cmake -S "$R/source/solvers/IncNavierStokesSolver/Tests/RANS" -B "$R/build" \
  -DCMAKE_BUILD_TYPE=Release -DRANS_BUILD_SOLVER=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  '-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math' \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_PREBUILT_SPATIAL_LIBRARY=/home/jay/build/wj-h14-mapped-backend-20261113/runtime/geometry/libSpatialDomains.so.5.9.0
cmake --build "$R/build" --target IncNavierStokesSolver \
  RANSWallResolutionProbe RANSWallMonitorChecks --parallel 8
```

Here R is the original H18 root. Three prequalification build versions and all
receipts remain. Final C++ source is 0209525d, independent of documentation HEAD.
Runtime setup clones the qualified private closure; do not substitute host libc.

`run_wj_h18.py` reserves every launched operation using the existing atomic
journal/process-group helpers. It admits only the frozen scientific matrix:

```sh
/usr/bin/python3 -B "$S/run_wj_h18.py" --root "$R" \
  --label target-A-v1 --seconds 750 --phase target --kind target \
  --attempts 4 --workers 2 -- \
  /usr/bin/python3 -B "$S/study_wj_h18.py" --root "$R" --mode target --case A
```

The receipt records the exact argv.
S is `$R/source/solvers/IncNavierStokesSolver/Tests/RANS`.
Target B similarly used 600 s/two attempts after measured admission, strictly
restoring A's unaudited absolute step2. The sole final strict operation used
`--phase target_readback --seconds 120` and `--mode strict`, after process exit.
No fresh target CP0 was created. Checkpoint suffixes are generation labels;
completion manifests and serialized absolute steps identify states.

The small sequence was legacy2, new-policy2, pending restart1, deliberate
post-implicit failure1. Its sealed dt=.01 remains a regression-fixture setting,
not a new target timestep. `wj_h18_checks.py` ran saved-small observer parity
and native guard/pending-integrity checks. The fourteen original selected Python
methods are enumerated in the focused-python receipt; three post-science
selector methods ran separately. The native protocol command was
`$R/build/RANSWallMonitorChecks $R/science/native-protocol-components`.
CTest registration alone is not execution evidence.

## Evidence and interpretation

Use the authoritative index, publication manifest and separately closed ledger.
Publication accounting is an OPEN snapshot, not the final closed journal.
The initial transport is a small exactly verified subset; local verification
lists remote-only originals. Original fields, load/row arrays, matrices and
histories stay Linux. Local summaries never replace producer/checkpoint pins.

The observer uses the sealed H16 final common support plus current witnesses.
It recomputes state fields, closure, current-K profiles/cutoff and reactions;
only immutable basis/geometry/reference data persist. It neither accumulates a
historical witness union nor rereads all archived ready packets. Existing final
ACK name validation remains; unrestricted archive/namespace cost is not qualified.

Online PASS, native scientific PASS and independent scientific PASS are separate.
Nondue audits are NOT_SCHEDULED. No all-time covariance/principal certificate,
unsaved extrema, exact integrated budget or nonlinear stress average follows.
Keep strong divergence, all 780 Q1 moments, frozen linear residual, genuine stage
defect and spatial residual distinct. The full-height inlet remains an idealized
control; physical 7 mm observations stay inside AWF. Long preservation, temporal
adequacy, physical height selection, terrain/hill runs and EP remain deferred.
