# CT-Accepted-6 evidence

The [report](../cao-tamura-accepted-6.md) distinguishes kernel/operator/static
qualification from the failed A95 evolving consistency gates. A93/E94E95
evolution transfers are conditional on complete A95 admission and have not
been attempted. Physical experimental applicability is a separate decision.

The [contract](contract.json) was committed before native execution. The
[pre-evolution addendum](pre-evolution-addendum-v1.json) also precedes every
attempted evolving step. [Execution notes](execution-notes-v2.json) retain the
diagnostic repairs. The [storage interruption](storage-interruption-v1.json)
records the user-coordinated mount recovery without deleting experiment data.
CT1--CT5 artifacts and the eight historical epsilon-order observations are
unchanged.

## Authoritative archive

All complete native packets, matrices, fields, histories, logs and timing data
remain on Linux:

```text
/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-6-20261113/recovered-v1
```

The original interrupted preparation remains in its parent directory. The
active filesystem is the user-mounted ntfs-3g volume. The old-path relocation
bridge and sealed CT4/CT5 runtimes have not been rewritten.

| Archive path relative to this root | Scope |
|---|---|
| `mathematics-v2/continuous` | Independent arbitrary-precision continuous fields and forcing, all models |
| `mathematics-v2/kernels` | Nine signed-curvature model cases and 15 boundary/domain controls |
| `native-operator-v1` | Small P3/P1 and P4/P2 all-row, weak, geometry and actual frozen-matrix checks |
| `initial-v1` | Four read-only projected initial representations; not reused as evolving histories |
| `targeted-v3`, `targeted-v4`, `static-independent-v2` | Zero-curvature/ground checks and compatible variable-height static collars |
| `evolution-v1` | Retained first flat-only ground-guard failure before any governing solve |
| `evolution-v2/A95-p4-dt002` | Successful actual trajectory, with retained original replay-diagnostic failure |
| `endpoint-supplement-v1`, `endpoint-classification-v2` | Read-only corrected replay of that unchanged trajectory |
| `evolution-v3` | Four remaining fresh A95 degree/timestep members through T=.016 |
| `physical-audits-v1` | Retained pre-repair scalar reconciliation attempt |
| `physical-audits-v2` | Native-only actual-stage diagnostics, explicitly insufficient for full quadrature qualification |
| `physical-audits-v3` | Complete four-control actual-state physical budgets |
| `primary-comparisons-v1`, `primary-comparisons-v2` | Original staged and complete-audit comparisons; failed gates stay failed |
| `controls-v1` | Retained ineffective observation-table orientation mutation |
| `controls-v2` | Corrected consumed-outlet-row negative controls |
| `regression-v1` | Fresh CT5 physical compatibility and curved restart/rollback/signature checks |
| `launch-*` | Closed owner-aware resource receipts, attribution and process CPU/elapsed costs |

## Commands and input identities

The exact commands, consumed configurations, executable and actual loaded-library
hashes are recorded in each `request-manifest.json`, execution receipt and
completed-artifact seal. Build requests identify the clean synchronized source;
`build-v12/results.json` records compiler settings and all frozen binaries.
The governing executable is byte-identical to `build-v8/runtime/IncNavierStokesSolver`;
later changes repair probe/reporting paths only.

These are the command forms used in the archive. Reproduction must use a new
isolated root and destinations. Do not run a writer against a sealed directory
or edit a checkpoint's method, time or signature.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export OMPI_MCA_orte_tmpdir_base=/tmp/nektar-ct6-ipc-20261113
export CT6ROOT='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-6-20261113/recovered-v1'
export CT6TEST="$CT6ROOT/source/solvers/IncNavierStokesSolver/Tests/RANS"
export CT6RUNTIME="$CT6ROOT/build-v12/runtime"
export LD_LIBRARY_PATH="$CT6RUNTIME/lib"
export TMPDIR="$CT6ROOT/tmp"

# Every numerical/build/analysis family uses one owned launcher.
python3 -B "$CT6TEST/ct_accepted3_command.py" --control "$NEW_CONTROL" -- \
  python3 -B "$CT6TEST/ct_accepted6_evolve.py" \
  --runtime "$CT6RUNTIME" \
  --contract "$CT6ROOT/source/docs/plans/cao-tamura-accepted-6-evidence/contract.json" \
  --admission "$CT6ROOT/admission-v5.json" \
  --model A95 --degree 4 --dt .001 --target .016 --output "$NEW_MEMBER"

# Complete read-only audit of an accepted native member.
python3 -B "$CT6TEST/ct_accepted3_command.py" --control "$NEW_AUDIT_CONTROL" -- \
  python3 -B "$CT6TEST/ct_accepted6_audit.py" \
  --runtime "$CT6RUNTIME" --member "$MEMBER" --output "$NEW_AUDIT"
```

The sequential `ct_accepted6_finish.py` driver reproduces the registered tests,
discriminating controls, flat/history controls and complete A95 endpoint audits.
Its root contains the completed primary trajectories with their original
identities. A collector exits 3 by design; a failed primitive admission exits 1.
Neither outcome may be interpreted as scientific admission merely because an
export completed. No physical state is changed by comparison or extrapolation.

Compact final classifications and archive hashes accompany the closed campaign
export. Large matrices and state arrays are deliberately excluded from Git.

## Completed decision

[qualified-v1/decision.json](qualified-v1/decision.json) contains the complete
eight-scope decision and raw primitive classifications. Its immutable artifact
seal and archive hash registry bind the full Linux packets.
[final-audit-v2/results.json](final-audit-v2/results.json) verifies every compact
sealed hash, recomputes admission and repeats the admission CLI. The CLI exits
1 for the retained degree/temporal failures. All five complete physical audits
and their 45 integration-uncertainty checks pass. A95/E94E95/A93 kernel and
static checks retain their declared scope; only A95 has evolving evidence.

The final verifier's own monitor closes after its result is sealed. Its
`resource-final.json` and `execution.json` are therefore copied alongside the
sealed files as separately identifiable closure receipts. The prior 25-family
inventory has no missing monitor closure; maxima are 14 project / 30 aggregate.
No continued numerical process remains from this campaign.

Additional exact command forms, each under a fresh owned control directory:

```bash
python3 -B "$CT6TEST/ct_accepted3_command.py" --control "$NEW_TEST_CONTROL" -- \
  ctest --test-dir "$CT6ROOT/build" -R '^CT6' --parallel 1 --output-on-failure

python3 -B "$CT6TEST/ct_accepted3_command.py" --control "$NEW_EXPORT_CONTROL" -- \
  python3 -B "$CT6TEST/collect_ct_accepted_6.py" \
  --root "$CT6ROOT" --output "$NEW_EXPORT"

python3 -B "$CT6TEST/ct_accepted3_command.py" --control "$NEW_GATE_CONTROL" -- \
  python3 -B "$CT6TEST/ct_accepted6_gate.py" \
  --input "$NEW_EXPORT/evolution-admission-input.json" --output "$NEW_GATE_RESULT"
```

The collector's expected exit is 3; the actual primitive gate's expected exit
is 1. The final verifier's zero exit means receipt verification passed, **not**
that the evolving formulation was admitted. The synthetic positive/negative
CLI tests and all 18 admission unit tests pass. The original ineffective
orientation fixture, interrupted/failed preparations, native-only receipts and
replay failures remain archived without retrospective threshold changes.
