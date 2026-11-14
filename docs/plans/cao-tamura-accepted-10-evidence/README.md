# CT10 compact evidence

`qualified-v1/decision.json` is the final primitive scope decision. It reports
passing refined-strip qualifications at .008/.02, the preserved baseline
pressure failure, passing target static checks, and failed required observer
capacity. Target trajectories are unattempted. Successful collection is not
physical admission.

## Records

- `contract.json`: prospective diagnostic, selection, history and physical gates.
- `single-contrast-v1.json`, `static-v1.json`: selected outer-band bisection,
  exact nodes, native counts, actual initial/boundary/reference hashes and
  zero-step admission, frozen before candidate trajectories.
- `target-preflight-v1.json`: prospectively bounded target stages/resources.
- `target-audit-capacity-v1.json`: explicitly non-evolving required-observer
  resource fixture, declared before its test.
- `qualified-v1/baseline-attribution-v3.json` and
  `candidate-attribution-v1.json`: complete signed component/band pressure
  attribution, native/independent reconstruction, pressure C/gauge/free-mass
  conditions, cancellation matrices, uncertainty and quadrature controls.
  Their `*-vectors.json` files retain the actual pressure coefficient vectors.
- `qualified-v1/*T008*.json`, `candidate-T020-v2.json`: complete derived field,
  traction, domain, mass, budget and timestep primitives with actual native
  report/history identities. `*-tractions.json` retains initial/endpoint physical
  profiles and all accepted-state traction summaries.
- `qualified-v1/common-space-v2.json`: equal-coordinate/weight/side comparisons,
  initial/evolved gradients, local constitutive inventories, regional pressure
  and reconstructed outlet stress. No coefficient-array spatial comparison.
- `qualified-v1/target-preflight-v1.json`: actual60,907-slot target initial,
  complete native-load, frozen-flow/action and quadrature results.
- `qualified-v1/target-audit-capacity.json`: actual initialized-field frozen-load
  capacity fixture; first native control fails `std::bad_alloc` at16GiB virtual
  address space. No actual consecutive-state physical budget is claimed.
- `qualified-v1/controls-v1.json`: exact face-only parity, original dense-export
  guard rejection and11 deliberately synthetic admission negatives. The signed
  load studies independently test omitted pressure/nonlinear-normal load and
  reversed natural traction.
- `qualified-v1/cost.json`, `resources.json`, `runtime.json`,
  `loaded-libraries.json`, `source-scope.json`: exclusive process timing,
  closed resource high-water records, exact executable/library/source identities.
- `qualified-v1/retained-attempts.json`: first diagnostic/controller preparations,
  with full artifact locations/hashes. No failed receipt is overwritten.
- `publication-v1/`: sealed post-collector provenance and its closed resource
  receipt. `publication-resource-v1/` is the publisher's separately closed outer
  receipt, copied after the publisher exited to avoid circular provenance.
- `observer-input-v2/`: separately sealed final observer-input supplement.
  Initial component vectors are identical to the original captures; baseline,
  refined-initial and refined-T=.02 loads agree with the old governing residual
  within 1.07e-14. The effective-input/lift policy is explicit, and all three
  tested raw/effective vectors coincide. No trajectory or earlier seal changed.
  Its build/test resources close and pass (14 project /18 aggregate peak).
  `observer-input-publication-resource-v1/` closes its outer publisher separately.

The compact Linux publication archive is:

```
/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-10-20261114/publication-v1.tar.gz
SHA256 8ded1fa7530114c10b6d0ec46889e25abc95a906dd3e83f0d1c45144bf679782
```

It contains only `qualified-v1` and `publication-v1`. The independently closed
publisher resource files are additional. Full fields, CSR matrices, checkpoints,
source bundles, native rows and rejected/preparation packets remain in the Linux
root. The target frozen field hash is
`60c53b83087a2c70b2a67faaf8cf09ef705ec39a04f04a5973f20244173368ff`.
Its capacity reproducer is `target-audit-capacity-v1/native/request.json` plus
`case.xml`, using the pinned CT9 probe. This input is a static capacity surrogate,
not a trajectory history.

The later read-only input supplement is archived separately:

```
/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-10-20261114/observer-input-v2.tar.gz
SHA256 f9eb327bb3b8ac87fe8731920f4d74e19a86e566b4fbc044bf32d653bfeefd77
```

It contains only the compact build/test results, immutable requests, seals and
their already closed resources. Its publisher receipt was copied after exit.
`files.json` binds all repository evidence, including the separate receipts.

## Reproduction

The exact argument arrays, inputs, loaded files and hashes are in each immutable
request/execution record. New output destinations are mandatory. Representative
commands, from a Linux shell using the tested source checkout, follow:

```bash
R='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-10-20261114'
C9='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-9-20261114'
C8='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-8-20261114'
S="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$R/tmp"
export OMPI_MCA_orte_tmpdir_base=/tmp/nektar-ct10-ipc

# Replay the exact original captured rows; no governing history is advanced.
python3 -B "$S/ct_accepted3_command.py" --control "$R/repro-load-resource" -- \
  timeout -k 30s 7200 python3 -B "$S/ct_accepted10_diagnose.py" \
  --ct8 "$C8" --ct9 "$C9" --runtime "$R/build-v2/runtime" \
  --contract "$R/source/docs/plans/cao-tamura-accepted-10-evidence/contract.json" \
  --capture "$R/diagnosis-v1" --output "$R/repro-loads"

# Reproduce the selected pair in fresh destinations using its sealed zero-step histories.
python3 -B "$S/ct_accepted3_command.py" --control "$R/repro-pair-resource" -- \
  timeout -k 30s 7200 python3 -B "$S/ct_accepted10_run.py" advance \
  --case candidate --target .008 --ct8 "$C8" --ct9 "$C9" \
  --runtime "$C9/build-v8/runtime" --preflight "$R/static-v1" \
  --contract "$R/source/docs/plans/cao-tamura-accepted-10-evidence/contract.json" \
  --output "$R/repro-pair"

# Minimal required-audit capacity failure; no target trajectory is generated.
python3 -B "$S/ct_accepted3_command.py" --control "$R/repro-capacity-resource" -- \
  timeout -k 30s 7200 python3 -B "$S/ct_accepted10_audit_capacity.py" \
  --ct9 "$C9" --runtime "$C9/build-v8/runtime" \
  --preflight "$R/target-preflight-v1" \
  --contract "$R/source/docs/plans/cao-tamura-accepted-10-evidence/target-audit-capacity-v1.json" \
  --output "$R/repro-capacity"
```

`ct_accepted10_run.py` returns nonzero on a failed physical pair.
`ct_accepted10_audit_capacity.py` returns nonzero on failed required capacity.
The collector may complete with an explicit failed/unattempted scientific branch;
`decision.json` retains the corresponding admission flag as false.

The final input-provenance control can be repeated under a fresh owned launcher
destination using `ct_accepted10_input_control.py --root "$R" --ct9 "$C9"
--runtime "$R/build-v4/runtime" --old-runtime "$C9/build-v8/runtime" --output
"$R/repro-input"`. Its compiled revision is `ec48cb35`; it advances no history.
