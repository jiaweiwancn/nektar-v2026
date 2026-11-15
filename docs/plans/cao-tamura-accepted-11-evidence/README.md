# CT11 evidence index

The observer passes old/new equivalence and target capacity. The actual target
pair reaches T=.002 and fails only `temporal:Pi`; longer continuation is also
outside the measured finite CPU/storage plan. See the [report](../cao-tamura-accepted-11.md)
and [sealed decision](qualified-v1/decision.json). Historical CT4–CT10 records
are unchanged. No terrain or experiment qualification is implied.

## Committed records

| Record | Purpose |
|---|---|
| [contract.json](contract.json) | Prospective scientific, observer, physical and resource gates |
| [resource-accounting-v2.json](resource-accounting-v2.json) | Prospective spent-plus-remaining closed-family accounting; v1 retained |
| [qualified-v1/](qualified-v1/) | Complete classified bodies, provenance, source scope, closed numerical resource receipts and exact request/seal |
| [summary-v1/results.json](summary-v1/results.json) | Small tables derived on Linux from the sealed collection and hash-checked timing traces |
| [native-receipts-v2/results.json](native-receipts-v2/results.json) | Direct extraction of retained native stage/trace/splitting fields; raw traction references kept separate |
| [pressure-failure-witness.json](pressure-failure-witness.json) | Paired pressure failure primitives and exact endpoint/history/stage bindings |
| [final-scope.json](final-scope.json) | Separate observer, physical pair, continuation-resource and terrain decisions |
| [commands.json](commands.json) | Exact completed owned-family commands and resource/source identities |
| [publication-v1/](publication-v1/) | Closed collector receipt, actual mount, compiler and dependency provenance |
| [resources-summary-v1/](resources-summary-v1/) and other resource directories | Separately closed reporting/extraction helpers; failed bookkeeping exits preserved |
| [files.json](files.json) | SHA256/size index of compact CT11 evidence; excludes itself |
| [publication-integrity.json](publication-integrity.json) | Strict JSON, immutable request/seal, lossless profile and source syntax checks; no new numerical conclusions |
| [next-package-prompt.md](next-package-prompt.md) | Target-pressure diagnosis and independent bounded throughput preflight |

All scientific artifact directories are immutable requests with completed file
seals. Verify every listed file SHA256 and the request SHA before using a result.
`qualified-v1` and `summary-v1` are separate sealed products. Final documentation,
this index and the failure witness supplement them; they do not edit their
classifiers. Some process exits are intentionally nonzero (target qualification,
negative native control, failed compile or bookkeeping). Worker compliance and
scientific statuses must be read independently.

Full coefficient vectors, checkpoints, native reports, matrices, stage packets
and raw audit outputs stay on Linux. Compact state hashes refer to the native
vector hashing scheme, while files have SHA256 identities. Do not treat them as
interchangeable. Complete actual mapped-library catalogs are keyed by digest in
`qualified-v1/loaded-libraries.json`; the recorded process mapping is stronger
evidence than an `ldd` query in a later collector environment.

The seven traction profiles in `qualified-v1/traction-profiles.json` are lossless,
deduplicated columnar tables. For each `ground`, `interface`, `interface_jump`
table, restore each original row with `dict(zip(table['columns'], row))` and
replace the table by that list. SHA256 of the restored profile's
`json.dumps(profile, sort_keys=True, separators=(',', ':')).encode()` must equal
its catalog key. All weights, coordinates, geometry/sides, gradients, endpoint
state, diffusivity, actual traction and row order remain present. No evolved
coefficient vector or full matrix is embedded in this catalog.

## Linux paths and identities

Verified mount: `/media/jay/新加卷1`, `/dev/sda2`, `fuseblk`.

```text
R=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-11-20261115
S=$R/source/solvers/IncNavierStokesSolver/Tests/RANS
CT10=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-10-20261114
CT9=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-9-20261114
GOV=$CT9/build-v8/runtime
PROBE=$R/build-v2/runtime
```

Starting revision: `8b2d56610ee12e8a2144a6becad2aa5150f263fb`, clean dedicated
Mac worktree, fetched remote matched. New probe built at `bfaa1617`, fresh target
driver at `3b2ba1c8`; publication at `eeed00f1`, tables at `047d7823`, retained
native receipt helper at `f11e5efe`. Full source revisions and file identities
are in each request. Governing executable SHA256 is
`5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7`;
probe SHA256 is
`20236b14eb4a90dea0f2a9a38b92e41135a83e9241c978fb265e185098f586e8`.
The governing runtime is unchanged, not a newly linked CT11 solver.

Immutable Linux locations:

| Directory | Scope |
|---|---|
| `$R/equivalence-v1` | 24 old/new controls; full all-row outputs external |
| `$R/capacity-v1` | Four target controls under 16-GiB RLIMIT_AS; static only |
| `$R/controls-v1` | 23 negative controls and native IBP/storage postchecks |
| `$R/target-v1/zero-0.001`, `zero-0.0005` | Independent genuine zero-step histories |
| `$R/target-v1/0.001`, `0.0005` | Actual accepted trajectories, native histories, stage packets, timing and four budgets |
| `$R/build-v1` | Failed warning-as-error build, retained |
| `$R/build-v2` | Passing isolated read-only probe build/runtime |
| `$R/qualified-v1`, `publication-v1`, `summary-v1`, `native-receipts-v2` | Sealed reporting/extraction products |
| `$R/resources-*` | Complete owned-family census, process CPU, wall, wait/suspension, PID/thread/resource data |

The original CT10 `target-audit-capacity-v1` failed allocation remains in its
original archive, bound through CT11's history/provenance. CT10's passing refined
strip and failed baseline pair are reused only within their unchanged scope.
CT9 operator/backend/restart/rollback controls are reused for the exact pinned
governing binary. No new native restart signature or numerical operator identity
was introduced by the observer.

## Reproduction

Use fresh destinations; completed requests cannot be overwritten. Activate the
authoritative Linux SDK, select the exact pinned governor/probe library path
for each process, then use the owned launcher. All exact arguments, timeouts,
resource roots and output destinations executed are in `commands.json` and the
underlying closed execution receipts. The final isolated source contains each
helper; rebuilding the probe at another revision requires fresh identities and
equivalence before physical use.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$R/tmp"
export OMPI_MCA_orte_tmpdir_base=/tmp/nektar-ct11-ipc
python3 -B "$S/ct_accepted3_command.py" --control "$R/resources-FRESH" -- \
  timeout -k 30s 7200 python3 -B "$S/ct_accepted11_capacity.py" \
  --ct9 "$CT9" --runtime "$PROBE" \
  --preflight "$CT10/target-preflight-v1" \
  --contract "$R/source/docs/plans/cao-tamura-accepted-11-evidence/contract.json" \
  --equivalence "$R/equivalence-v1" \
  --output "$R/capacity-FRESH"
```

The example reruns the exact capacity fixture into a fresh destination. Target stepping additionally
requires the complete positive equivalence, capacity, verification and closed
resource receipts, exact CT10 preflight/strip requests, contract and pinned
runtime. A successful observer command alone never admits evolution.

Publication archive hashes (archives remain outside Git):

| Linux file | SHA256 |
|---|---|
| `$R/publication-v1.tar.gz` | `7819a0f1befe09a1a130acff682ae75d3b1247bffb516f57524e8a75fb6810c7` |
| `$R/summary-package-v1.tar.gz` | `5aa755d8e1a257aff193685c2803d3bada79c3b5eb04a27abc7a047c2b6639a1` |
| `$R/native-receipts-package-v1.tar.gz` | `38869f5410bb680332857e5e1f91e554bc8a3e956dfc2916ac0427aad8df742e` |

Publisher/summary/package resource receipts close after each archive is written
and are committed separately, avoiding a circular self-seal. A bookkeeping
bundle-ref error and failed helper launch precede the successful
`native-receipts-v2`; no numerical state was created by those failures.
