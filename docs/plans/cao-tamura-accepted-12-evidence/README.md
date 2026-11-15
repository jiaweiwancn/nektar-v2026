# CT12 evidence and reproduction

[Report](../cao-tamura-accepted-12.md), [contract](contract.json),
[method](method.md), [evaluated decisions](evaluated-v1/decision.json),
[final closed accounting](closed-publication-accounting.json),
[tested source scope](source-scope.json),
[next prompt](next-package-prompt.md).

The Linux archive is
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-12-20261115`.
The verified mount is `/dev/sda2`, fuseblk. Large native packets, matrices,
pressure vectors, fields, R3 histories and raw logs stay there. No CT4–CT11
archive was pruned or rewritten. This export contains complete compact result
bodies, primitive classifications, selected witnesses, source/runtime identities
and hashes; a successful collector exit does not mean physical admission.

## Authoritative inputs and runtimes

- CT8 physical manifest SHA256:
  `004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296`.
- Governing runtime: `ct-accepted-9-20261114/build-v8/runtime`, executable SHA256
  `5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7`.
- Probe runtime: `ct-accepted-11-20261115/build-v2/runtime`, executable SHA256
  `20236b14eb4a90dea0f2a9a38b92e41135a83e9241c978fb265e185098f586e8`.
- Component: `CSR1:EigenLU-LDLT-IR3-F40M:FGMRES80x1000`.
- The pinned MPI-enabled builds execute one serial rank; numerical libraries
  use one thread. No CT12 C++ rebuild or shared-install change occurred.
- Actual loaded native libraries are in the replay/probe execution receipts;
  actual Python/SciPy mappings are in the axis receipt. Solver and diagnostic
  source identities are distinct. First replay-step temporary RSS/map telemetry
  was lost before a wrapper receipt; second-step telemetry and all exact stage
  data survived. This limitation is not replaced by an invented measurement.

## Newly generated and reused data

| Archive subdirectory | Scope |
|---|---|
| `pressure-v4` | Eight native aggregate state/stage exports plus held q48/q64 controls; q64 fails allocation under16GiB |
| `pressure-v5` | Completed six-correction algebra and signed attribution; exact hash-checked reuse of those packets; no new native export |
| `profile-v1/step1`, `step2`, `rollback` | Exactly two accepted replay operations and one deliberate rejected operation, with original runtime and genuine histories |
| `profile-v3` | Completed final four-control physical audit and lossless-retention proof; no native operation repeated during wrapper recovery |
| `retention-v1` | Content-addressed lossless blobs, required-evidence proof and exact restoration index |
| `axes-v2` | Completed bounded-factor pressure/axis/ground/interface/feedback diagnosis; its inner external-traction postprocessor is explicitly invalidated |
| `axes-v3` | Correct radial external traction; retained valid axis results; independent old-implementation traction comparison |
| `face-v1` | Independent physical inlet/outlet shear Stokes identity and pressure-image/orientation test |
| `static-v1` | One fresh zero-time first-inner-slab contrast;62533 slots, blocked by40497972 predictor factor entries |
| `evaluated-v1` | Primitive recomputation, negative controls, resource/campaign decisions and compact bodies |

Retained CT11 target reports, initial/accepted/flow/pre-k states, CSR blocks,
free/prescribed IDs, histories and their hashes are reused. CT10 target static
rows/inventories are reused on the identical original mesh. All native
coordinates, weights, sides and generating owners remain those of the pinned
runtime. Interface flux observations are not duplicate governing loads.

The old CT11 physical result remains failed only at `temporal:Pi`. No new target
pair, timestep halving, evolved contrast, longer continuation or terrain run was
performed. The actual physical horizon remainsT=.002.

## Preserved failures and correction notices

- `pressure-v1`: tiny-RHS relative-only diagnostic work bound; owned interrupt.
- `pressure-v2/v3`: capped small-RHS diagnostic stopping flags, followed by
  prospective primitive/tighter-repeat admission amendments. They are not
  silently converted to completed experiments.
- `pressure-v4`: late q64 native aggregate `std::bad_alloc`; reconstruction data
  recovered without repeating native exports.
- `profile-v1/v2`: checkpoint-directory and empty rejected-checkpoint wrapper
  errors after completed native work; fixed without repeating accepted steps
  or the rejection.
- `axes-v1`: successive diagnostic factor reference cycles hit16GiB; explicit
  factor release completes the same reconstruction within the original bound.
- `axes-v2`: ground D(0) was wrongly used at finite-height inner external faces.
  [The correction](external-traction-correction.json) invalidates only those
  postprocessed values and binds the corrected receipt to its sealed parent.
- `static-v1`: the original40M fill guard stops the selected mesh. Missing later
  controls are unattempted, not failed mathematical reference tests.

## Commands

Every numerical/control command is recorded exactly in
`evaluated-v1/resources.json`, including closed CPU, elapsed time, suspension
information, PIDs, thread counts and high-water attribution. The final
`prepublication-accounting.json` adds later closed inspection/collection families,
and `publication-execution.json` closes the final packaging family. Their additive
`closed-publication-accounting.json` records47 closed families,9598.852857 owned
process-CPU seconds and7650.416142 summed active-wall upper seconds. The package
worker peak is14 project/18 aggregate. No numerical family remains open.
Keep source/Git
work on Mac; execute the following only in a fresh isolated Linux destination.
Do not reuse a completed destination label.

```bash
export CT12_ARCHIVE='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-12-20261115'
export CT_RUNS='/media/jay/新加卷1/nektar-v2026-runs'
export CT_TESTS="$CT12_ARCHIVE/source/solvers/IncNavierStokesSolver/Tests/RANS"
export CT_FRESH="$CT_RUNS/ct12-reproduce-NEW"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$CT_FRESH/tmp"
export OMPI_MCA_orte_tmpdir_base=/tmp/nektar-ct12-ipc
mkdir -p "$TMPDIR"
```

Recompute the read-only diagnosis from the retained identical native packets:

```bash
python3 -B "$CT_TESTS/ct_accepted3_command.py" \
  --control "$CT_FRESH/resources-pressure" -- \
  timeout -k 30s 2600 python3 -B "$CT_TESTS/ct_accepted12_pressure.py" \
  --ct10 "$CT_RUNS/ct-accepted-10-20261114" \
  --ct11 "$CT_RUNS/ct-accepted-11-20261115" \
  --runtime "$CT_RUNS/ct-accepted-11-20261115/build-v2/runtime" \
  --contract "$CT12_ARCHIVE/source/docs/plans/cao-tamura-accepted-12-evidence/contract.json" \
  --reuse-packets "$CT12_ARCHIVE/pressure-v4" --output "$CT_FRESH/pressure"
```

A separate reproduction can execute the bounded replay using
`ct_accepted12_profile.py --ct10 ... --ct11 ... --runtime <CT9/build-v8/runtime>
--probe <CT11/build-v2/runtime> --contract ... --output <fresh/profile>` under the
same owned launcher. Its default is two accepted steps plus one declared
rejection. The recorded `--recover/--finish-only` commands reproduce the actual
controller recovery lineage; they never repeat existing native work.

Use `ct_accepted12_axes.py` for axis/traction reconstruction,
`ct_accepted12_face.py` for the independent face identity and
`ct_accepted12_static.py` with `static-first-inner-v1.json` for the one selected
zero-time contrast. The static driver deliberately returns a collection result
after a native guard; inspect its primitive classification, not its exit code.
The final `ct_accepted12_collect.py` command binds pressure, corrected axes,
profile, face and static paths and recomputes admission.

[Publication hashes](publication-provenance.json) and the compact-file checksum
[index](compact-checksums.json) bind the exported bodies to exact archive files.
The transfer archive SHA256 is in `transfer-sha256.json`. The evaluated seal is
copied byte-for-byte; publication metadata is additive. The stdlib-only
`publication-driver.py` is the exact final Linux packaging source, with its
command and closed resource receipt retained. Mac performs transport verification
and addition of those closed receipts, without producing new numerical claims.
Full pressure vectors,
all native requests, actual stage packets and failure logs remain at the listed
archive paths for independent reconstruction.
