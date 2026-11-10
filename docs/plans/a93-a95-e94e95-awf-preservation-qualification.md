# Native audit and neutral-reference preservation qualification

Starting revision: `93b6c6b31a58a976ce4462526114243646e6552d`, clean local branch
`codex/a93-a95-e94e95-awf`; remote fetch agrees. Previous primary native binary:
`1950dedf4e04a903a3e919b64bc3ce7a25c25358`. The numerical source hashes at the
reviewed head agree with that binary's manifest. Its immutable private runtime
is reusable for before/after controls, not relabelled as a new build.

Mac source worktree remains `.worktrees/a93-a95-e94e95-awf`. Linux work is under
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-preservation-20261110`, with
separate source, build, runtime and experiment subdirectories. The older A93
and active CLS directories are read-only dependencies/evidence. Initial Linux
conservative process census: 18 relevant threads. Build allocation: eight;
native jobs initially serial, MPI ranks 1, all numerical-library thread counts
1. Additional concurrency requires a fresh census and remains below 72.

## Contract and defects

The first-work contract and supplied TeX specifications govern unchanged.
Their SHA-256 remain respectively
`056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64` and
`17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40`.
No continuum coefficient, source, wall calibration, epsilon semantics or
one-sweep ordering is changed. No anchoring, cancellation forcing, field reset,
radicand clipping or new positivity device is permitted.

Audit corrections:

* Required cutoff checks distinguish accepted budget profiles, active EP
  corrections and inactive EP baselines. Actual scalar/constitutive/covariance
  checks remain independent. Activity, fallback and limiting reason are logged.
* New-model forward audits omit unsupported tangent stencils. The exact A95
  near-boundary physical witness goes through the production diagnostic helper;
  unavailable tangent evidence is explicitly not an ellipticity certificate.
* All covariance callers use one scale-normalized symmetric eigensolver,
  including an exactly zero branch. Hard covariance eligibility uses the same
  relative reference-k tolerance under coherent rescaling.

## Predeclared numerical gates

Each result distinguishes execution, expected rejection, analytical identity,
representation accuracy, quadrature uncertainty, spatial convergence and
temporal order. Missing coverage or a failed required check returns a nonzero
process status after available reports and reproducing XMLs are preserved.

Analytical weak balances use native test-function values/gradients and physical
integration points, with independently evaluated physical constitutive laws and
wall profiles. Scaled weak defects must be below 1e-9 (absolute unscaled floor
1e-12). Shared interface rows occur once. Accepted-wall source and diffusion
moments are retained separately; they need not separately vanish.

Frozen-state quadrature series: wall-normal 32/48/64, tangential and outer-volume
8/12/16, corresponding mass and interface 8/12/16. The coefficient state and
polynomial spaces stay identical. Radial reference projection uses its own
reported native Gauss32 rule and geometric subpanels. Require the last
quadrature change below 1% of an interpreted spatial-error change, with a
scaled absolute floor 1e-10; otherwise classify the spatial conclusion as
quadrature unresolved. Common evaluation-grid uncertainty is checked separately.

Independent h series: three nested radial levels at p=3, with declared physical
breakpoints (cutoff/intercept/roughness/H) aligned in the base mesh. Independent
p series: p=3,4,5 on one fixed mesh, with pressure degree p-2. Values, gradients,
alpha=0 spatial residual groups and pointwise interface quantities use fixed
physical scales. A resolved observable requires sustained reduction at both
increments (ratio <=0.85), beyond five times quadrature/evaluation uncertainty.
Near-zero observables are reported as floors, not convergence orders.

Temporal triplets use dt=1e-4,5e-5,2.5e-5 at T=.0016 with identical initial
coefficient states, fixed quadrature/spaces and consistent native histories.
Comparisons use matching physical grids and weighted RMS, plus one common
global pressure gauge. Differences above five times their uncertainty use the
initial first-order band 0.8<=p<=1.2; smaller signals are unresolved. A targeted
finer dt or bounded extension may diagnose a specific failed/unresolved order.
Perturbation response is measured against the unperturbed numerical trajectory
at the same time, separately from reference drift; no monotone-decay criterion
is imposed.

The six natural-kappa model/wall combinations receive temporal triplets. All
twelve model/wall/policy configurations receive selected spatial and short
native checks; prescribed-kappa temporal coverage is explicitly tabulated.
Capabilities not exercised are not inferred from successful endpoints.

## Results

Completed through bounded native consistency qualification. The final machine receipt is unresolved only for the below-floor unperturbed epsilon orders; see the numerical tables and evidence index below.

## Implementation and source mapping

| Requirement | Authoritative implementation / verification |
|---|---|
| Formulation/activity-aware cutoff eligibility | `RANSWallEvolutionFeasible.cpp`: separate informational and required budget minima, active/inactive/fallback counts; `RANSEquilibriumWallEvolution.cpp`: specific rejection reasons |
| Forward validity versus unavailable tangents | `RANSWallForwardDiagnostics.h`, used by `RANSWallEvolutionDiagnostics.cpp`; exact witness in `RANSNativeDiagnosticChecks.cpp`, full native polynomial witness in `validate_rans_native_audits.py` |
| Scale-independent covariance diagnostics | One normalized eigensystem in `RANSWallSolveDiagnostics.h`; all minimum-eigenvalue/eigenvector callers use it; the original tensor is unchanged |
| Actual A95 domain margin | Forward constitutive output in `RANSModels.*`; native inventory/range and limiting-point receipts in `RANSWallEvolutionFeasible.cpp` |
| Strict execution, coverage, sample identities and pressure gauge | `rans_qualification.py`, matrix runner, native qualification and collectors; declared required check names cannot be removed to make coverage pass |
| Analytical weak oracle | `rans_preservation_oracle.py`: independent dimensional stresses, sources and wall primitives, evaluated against native exported test functions and boundary traces |
| Unanchored alpha=0 replay / common physical grid | `RANSPreservationProbe.cpp`, `RANSWallEvolutionPatch::PreservationPacket`; descriptor, mesh/space and boundary identity checks; existing legacy replay guard retained |
| Fixed-coefficient quadrature histories | Optional **fresh time-zero initialization artifact** from the probe and `check_rans_frozen_initial_history.py`; native IMEX initialization and native checkpoint writer/signature, no governing advance and no evolved-history conversion |
| Quadrature, h, p and temporal jobs | `qualify_rans_preservation.py`, `check_rans_preservation_uncertainty.py`, `reanalyse_rans_temporal.py`; execution, accuracy, convergence and expected rejection are separate receipts |
| Small CTest controls | Native audit fixtures, normalized spectrum/domain checks and Python failure-contract tests; larger campaigns remain explicit jobs |

The new history-export option is restricted to admissible identical state and
reference coefficient bytes at time zero, native IMEX1 and the qualified
unforced neutral boundary contract. It constructs fresh native solution and
explicit-derivative history slots and the **target** session's full physics
signature. The ordinary solver checks that signature when loading it. A changed
descriptor/identity, nonzero time or unequal state/reference is rejected. No
signature in an existing checkpoint is edited. Both default replay and export
leave the supplied physical coefficients unchanged and perform zero steps.

## Observed discrepancies and preserved investigations

1. The old inactive EP hard gate combined unused corrected-profile budgets with
   actual eligibility. Positive scalar/covariance baselines were consequently
   rejected. Required-cutoff eligibility now follows the actual formulation and
   correction activity. Scalar and constitutive/covariance requirements remain
   mandatory.
2. The old diagnostic stencil left A95's radicand domain at a physically valid
   state. The exact dimensional witness has radicand about `5.68e-15` in the
   production long-double calculation, positive covariance, and two invalid old
   stencil probes. A trace-compatible native polynomial witness independently
   reaches a minimum physical margin about `4.23e-14`; the old executable rejects
   it, while the corrected native audit accepts it. Tangent availability is zero,
   not an ellipticity certificate.
3. The old dimensional eigensolver misses the off-diagonal entries of the
   specified tiny tensor. Its reported minimum is `+6.67e-17`, instead of
   `-3.33e-17`. A separate full native experiment coherently scales k, epsilon,
   velocity, reference friction velocity and molecular viscosity: the old hard
   audit rejects the ordinary-scale indefinite state but accepts its `k=1e-16`
   counterpart. The corrected hard audit rejects both specifically for
   covariance. Actual scalar and cutoff eligibility are valid in that fixture.
4. The first strict temporal reporting attempt omitted weight identity on the
   analytical comparison object. It correctly failed, although its native
   trajectories completed. Original reports and failed driver receipts remain
   under `preservation-temporal-Q0-*`. Read-only reanalysis validates final
   coefficients, initial reference, physical coordinates, weights and terminal
   time before recomputing metrics. It does not rerun or alter those trajectories.
5. Accepted `ReferenceK(r)` itself integrates with the wall normal rule, before
   the fixed Gauss32 radial projection. Thus independently regenerating the
   initializer at normal orders 32 and 64 changes some coefficient values by up
   to `1.75e-13`. Those original endpoint controls fail the exact initial-byte
   comparison and are retained in `uncertainty/`. The strict frozen residual
   study already reused identical coefficients. Fresh, validated zero-step
   histories additionally remove this initialization difference from the
   endpoint quadrature controls. The first exporter correctly triggered the
   native signature guard because it omitted the template's `Chart` function
   convention; its receipts are retained, and the exporter now uses the shared
   native signature with the actual convention.
6. A 0.5% initial epsilon perturbation resolves five of the six natural-kappa
   epsilon signals. A95/EP remains just below the declared floor. One targeted
   1% perturbation triplet resolves that signal without changing the closure,
   source, boundary data, time horizon or acceptance band. Both amplitudes and
   the unresolved original control are retained.

## Interpretation of the measurements

Analytical fields and analytical gradients are distinct from the constrained
polynomial reference and the evolved numerical state. The independent weak
oracle includes exterior fluxes, outer epsilon boundary fluxes, complete normal
and shear traction, both pressure regions with one gauge, and the shared k row
once. Accepted-wall source moments are appreciably nonzero; their diffusion
moments cancel them. EP plateau source moments are zero pointwise. No independent
inner epsilon equation or extra two-sided epsilon-flux condition is introduced.

Nektar's flat 2D `(x,y)` session coordinates correspond to the notes' `(x,z)`;
exported physical tensors and coordinates retain all three components. The
common evaluation grid uses 204 fixed radial panels, GL10 in the normal direction
and GL8 horizontally. It contains every compared element boundary and all
initial reference breakpoints. Native meshes stay fixed during advancement;
profile quadrature panels may follow the current cutoff. Cutoff ranges and
all accepted/stage extrema are retained separately from final-state extrema.

The global flat coefficient mass uses native Gauss `max(8,p+3)`, independently
of the optional patch diagnostic/regional override. For p=3,4,5 its integrands
have degree at most ten and Gauss8 is already exact in exact arithmetic.
Outer volume, source/profile, interface and diagnostic/regional rules are
explicitly recorded. The accepted reference's nested integral order is also
recorded; it is not hidden behind the Gauss32 projection label.

All comparisons use weighted RMS, `sqrt(integral(error^2)/volume)`, and separately
report true L2. A single volume-weighted global Pi gauge is removed, preserving
the pressure jump between regions. Array lengths, coordinates, weights,
regions, field identities and terminal times are checked before comparison.
GL10 also integrates squared differences between the compared p<=5 polynomial
states exactly on the common element partition in exact arithmetic. GL14/GL10
re-evaluation independently checks the nonpolynomial analytical-reference error.

Spatial residuals are native **alpha=0** residuals with unchanged physical
scales. Their momentum components, global/inner/outer/shared-interface k,
outer epsilon, continuity, gauge, dependent trace and pointwise traction/flux
receipts are separate. The evolved-state alpha=0 receipts are distinct from
small time-discrete stage residuals containing mass/history terms. A finite
polynomial approximation to the logarithm/piecewise k is not a roundoff
stationary state; the substantial coarse residuals are reported, not canceled.

Temporal orders use identical initial coefficient hashes within each triplet
and T=0.0016. Every order above five times the declared `1e-10` uncertainty floor
must lie in [0.8,1.2]. Below-floor values remain unresolved. Reanalysis returns
exit **2** for that outcome, exit **1** for failed/missing required checks and
exit **0** only for resolved passing gates. The refreshed epsilon trace lag is
reported as splitting error and scales approximately with dt; no evolved trace
is repaired after the k solve. Perturbation response is measured against the
unperturbed numerical solution at the identical time, separately from drift
relative to the analytical reference; no monotone-decay assertion is imposed.

## Numerical results

Final required numerical gates pass; the overall receipt is **unresolved / exit 2**, with no failed required checks, because eight unperturbed outer-epsilon temporal signals are below the predeclared resolution threshold. Their numerical orders remain null. All six natural-kappa model/wall combinations additionally have resolved perturbed epsilon orders.

| Evidence | Result |
|---|---|
| Native audit defects | Nine before/after controls; all intended outcomes observed |
| Spectrum / domain helper | 696 checks; normalized spectra, eigenvectors, hard decisions and exact A95 witness pass |
| Independent model/wall kernels | 1346 checks pass |
| Native weak oracle | 14 evaluations, including all twelve configurations; worst scaled required defect 3.11e-13 |
| Spatial study | 34 initialized/frozen checks; independent h and p gates pass |
| Temporal study | Twelve unperturbed triplets (36 runs); 52 resolved field orders pass, eight epsilon orders unresolved |
| Epsilon signal controls | Six natural-kappa combinations pass; 0.5%, with one targeted 1% A95/EP triplet |
| Native matrix | Twelve combinations, 16 steps each; previous binary checkpoint states/histories bitwise identical |
| Regression | 22 continuation/guard controls, A95 cache, reviewed CLS restart and 15 distinct CTest targets pass |

Frozen residual changes from the final two quadrature rules are below 1e-13 for the principal groups and below the declared uncertainty criterion. The largest measured Q2/Q0 endpoint field RMS difference is 7.41e-11 (pressure); the native starting coefficients are byte-identical. Common-grid GL10/GL14 changes in analytical value/gradient error estimates are below 9e-17. No bound projection activates in these trajectories.

### Separate h and p results

E94/E95, natural kappa. h0/h1/h2 use p=3; p4/p5 use the fixed h0 mesh. Momentum and scalar residuals are dimensionless mass-dual alpha=0 norms on unchanged physical scales.

| Wall | Space | U value RMS | U gradient RMS | Momentum x | Inner k | Outer epsilon |
|---|---|---:|---:|---:|---:|---:|
| EP | h0 | 0.0003255 | 0.03334 | 3.907 | 0.02257 | 0.0017 |
| EP | h1 | 2.371e-05 | 0.006378 | 1.183 | 0.00896 | 0.000409 |
| EP | h2 | 1.689e-06 | 0.0009425 | 0.2512 | 0.001017 | 9.25e-05 |
| EP | p4 | 5.238e-05 | 0.008023 | 0.8063 | 0.01573 | 0.0001692 |
| EP | p5 | 6.531e-06 | 0.001612 | 0.1081 | 0.001246 | 1.681e-05 |
| Accepted | h0 | 0.0003305 | 0.03383 | 4.042 | 2.796 | 0.0017 |
| Accepted | h1 | 2.408e-05 | 0.006426 | 1.208 | 1.13 | 0.000409 |
| Accepted | h2 | 1.715e-06 | 0.0009534 | 0.2565 | 0.2491 | 9.25e-05 |
| Accepted | p4 | 5.319e-05 | 0.007976 | 0.7728 | 0.6386 | 0.0001692 |
| Accepted | p5 | 6.632e-06 | 0.001636 | 0.1106 | 0.1008 | 1.681e-05 |

Accepted k value RMS decreases `3.12e-5 -> 4.05e-6 -> 3.32e-7` under h refinement; its gradient RMS decreases `1.55e-2 -> 2.91e-3 -> 4.29e-4`. At fixed mesh p=3,4,5 the corresponding value errors are `3.12e-5, 8.30e-6, 1.41e-6` and gradient errors `1.55e-2, 3.62e-3, 7.84e-4`. EP constant-k representation is at roundoff, while its native k residual still reflects the velocity-gradient representation error.

| Wall | Shared k row h0 -> h1 -> h2 | Outer k h0 -> h1 -> h2 |
|---|---:|---:|
| EP | 0.0002623 -> 1.7e-05 -> 1.001e-06 | 1.617e-05 -> 1.308e-06 -> 1.186e-07 |
| Accepted | 0.0002906 -> 2.141e-05 -> 1.339e-06 | 1.617e-05 -> 1.308e-06 -> 1.186e-07 |

A93/A95 selected normal-stress checks (h0 -> h1, fixed p=3):

| Model / wall | Vertical momentum residual | Full pointwise traction-jump RMS |
|---|---:|---:|
| A95 / EP | 0.0002655 -> 6.294e-05 | 4.401e-14 -> 9.131e-14 |
| A95 / Accepted | 0.0002655 -> 6.294e-05 | 4.752e-14 -> 9.643e-14 |
| A93 / EP | 0.0006755 -> 0.0001602 | 3.829e-14 -> 7.842e-14 |
| A93 / Accepted | 0.0006755 -> 0.0001602 | 3.761e-14 -> 7.811e-14 |

The pointwise traction jumps are already below the absolute floor; they are not claimed to converge from roundoff. Their small size does not replace the separate nonzero shared weak-row residual. All twelve model/wall/policy configurations have the selected coarse/fine spatial receipts in `spatial-results.json`.

### Temporal orders at T=0.0016

`dt=1e-4, 5e-5, 2.5e-5`; identical initial coefficient hashes within each triplet. `w` is the native session variable `v`, i.e. physical vertical velocity.

| Model | Kappa | Wall | p(U) | p(w) | p(k) | p(epsilon) | p(Pi) |
|---|---|---|---:|---:|---:|---:|---:|
| E94E95 | Natural | EP | 0.978 | 0.978 | 0.968 | unresolved | 0.992 |
| E94E95 | Natural | Accepted | 0.978 | 0.978 | 0.974 | 0.996 | 0.992 |
| E94E95 | Prescribed | EP | 0.978 | 0.978 | 0.968 | unresolved | 0.992 |
| E94E95 | Prescribed | Accepted | 0.978 | 0.979 | 0.974 | 0.996 | 0.992 |
| A95 | Natural | EP | 0.992 | 0.963 | 0.964 | unresolved | 0.995 |
| A95 | Natural | Accepted | 0.978 | 0.962 | 0.974 | unresolved | 0.994 |
| A95 | Prescribed | EP | 0.987 | 0.963 | 0.970 | unresolved | 0.994 |
| A95 | Prescribed | Accepted | 0.978 | 0.962 | 0.974 | unresolved | 0.994 |
| A93 | Natural | EP | 0.978 | 0.945 | 0.967 | unresolved | 0.993 |
| A93 | Natural | Accepted | 0.978 | 0.945 | 0.974 | 0.991 | 0.993 |
| A93 | Prescribed | EP | 0.978 | 0.946 | 0.968 | unresolved | 0.993 |
| A93 | Prescribed | Accepted | 0.978 | 0.946 | 0.974 | 0.992 | 0.993 |

Resolved unperturbed orders range from 0.945 to 0.996. The unresolved entries are all six EP epsilon cases and both A95 accepted-wall epsilon cases. No prescribed-kappa temporal coverage is inferred: all six prescribed-kappa configurations actually have their own triplets.

| Natural-kappa signal control | Initial epsilon amplitude / epsilon(H) | Measured epsilon order |
|---|---:|---:|
| A93 / Accepted | 0.5% | 0.995 |
| A93 / EP | 0.5% | 0.993 |
| A95 / Accepted | 0.5% | 0.999 |
| A95 / EP | 1.0% | 0.998 |
| E94E95 / Accepted | 0.5% | 0.998 |
| E94E95 / EP | 0.5% | 0.997 |

The largest refreshed trace lag decreases roughly by four when dt decreases by four: for E94/E95 it is `9.62e-8 -> 2.44e-8` (EP) and `1.35e-7 -> 3.41e-8` (accepted). A93/A95 show the same scaling. The raw lag, source-splitting receipts, final alpha=0 spatial residuals and time-discrete stage residuals are all retained. These are first-order consistency measurements on the fixed qualified spaces, not a complete temporal-order or long-time stability campaign.

The initial accepted-wall coefficient differences observed in the uncontrolled quadrature comparison are `4.77e-14` (E94/E95), `1.19e-13` (A95) and `1.75e-13` (A93). The final endpoint controls eliminate them using fresh native histories. Their native physics signatures and starting-state byte hashes are verified; all four new replay/export rejection controls pass.

Serial runs through the MPI-enabled SDK emit the same PMIX dstore permission warning seen in the earlier baseline. Exit statuses, completed histories and numerical checks are recorded independently. No MPI numerical qualification is inferred from these runs.


The compact receipt index is [preservation-evidence/README.md](a93-a95-e94e95-awf-preservation-evidence/README.md); [overall-gates.json](a93-a95-e94e95-awf-preservation-evidence/overall-gates.json) gives the machine decision.

## Provenance and reproduction

The primary frozen solver was compiled from numerical C++ at `3655ac3` and frozen
with source manifest `0db2f1f`. Subsequent production numerical C++ is unchanged.
Later verification-only executables have their own manifests. The primary solver SHA-256 is
`cb08d0c4689177913960932f3e6a6c3cca6fc073937c8b6f5240137541d72e0a`.
The exact source files, compiler/options, executable and linked-library hashes, actually loaded
libraries/plugins, configurations and thread environment accompany the compact
receipts. Direct linked dependencies are private copies; dynamically loaded
system MPI plugins are identified by path and hash. This is not a claim of a
containerized host environment.

The original and new runs, including failed inputs, checkpoints, field files,
geometry/basis packets and timing logs, remain outside the checkout under the
Linux experiment root. Native states/histories in all twelve original short
configurations match the previous `1950dedf` binary's checkpoint payloads
bit-for-bit. Current restart/cache/rollback and standard/CLS controls are rerun;
legacy CLS continuation from the reviewed baseline is also compared through
binary checkpoint payloads rather than JSON representations of signed zero.

Representative commands (Linux, empty destinations for new runs; existing
receipts are never overwritten):

```bash
export RANS_PRESERVE_ROOT=/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-preservation-20261110
export RANS_TESTS="$RANS_PRESERVE_ROOT/source/solvers/IncNavierStokesSolver/Tests/RANS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$RANS_PRESERVE_ROOT/tmp"
export LD_LIBRARY_PATH="$RANS_PRESERVE_ROOT/runtime/lib"

cmake -S "$RANS_TESTS" -B "$RANS_PRESERVE_ROOT/build" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_BUILD_SOLVER=ON \
  -DRANS_PREBUILT_SPATIAL_LIBRARY="$RANS_PRESERVE_ROOT/dependency/libSpatialDomains.so.5.9.0" \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$RANS_PRESERVE_ROOT/build" --parallel 8 --target \
  IncNavierStokesSolver RANSPreservationProbe RANSNewModelsProbe \
  RANSNativeDiagnosticChecks RANSKernelProbe RANSConstitutiveProbe \
  RANSAcceptedAWFProbe RANSProfileBranchChecks

python3 -B "$RANS_TESTS/validate_rans_native_audits.py" \
  --solver "$RANS_PRESERVE_ROOT/runtime/IncNavierStokesSolver" \
  --descriptor-probe "$RANS_PRESERVE_ROOT/runtime/RANSNewModelsProbe" \
  --output NEW_AUDIT_DIRECTORY
python3 -B "$RANS_TESTS/qualify_rans_preservation.py" \
  --root "$RANS_PRESERVE_ROOT" --phase spatial
for model in E94E95 A95 A93; do
  python3 -B "$RANS_TESTS/qualify_rans_preservation.py" \
    --root "$RANS_PRESERVE_ROOT" --phase temporal --models "$model" \
    --label "Q0-$model" --quadrature-index 0
  # Exit 2 records unresolved signals; inspect the receipt and continue other jobs.
done
python3 -B "$RANS_TESTS/check_rans_preservation_uncertainty.py" --root "$RANS_PRESERVE_ROOT"
python3 -B "$RANS_TESTS/check_rans_frozen_initial_history.py" --root "$RANS_PRESERVE_ROOT"
python3 -B "$RANS_TESTS/collect_rans_preservation.py" \
  --root "$RANS_PRESERVE_ROOT" --output NEW_COLLECTION_DIRECTORY
```

These commands are the recorded experiment's layout, not an instruction to
replace a shared installation. To repeat the trajectories, create a new source
worktree/build/runtime/run root first. The private geometry dependency and frozen
runtime identities are in the manifests. Freeze the built executable and direct
libraries with the repository's `continue_cls_epawf_pimple.freeze` helper before
running. Larger campaign directories deliberately reject reuse; use their
preserved XML/request files or a new root. Python >=3.10, NumPy and SciPy are
required explicitly at configuration and command-line entry. Fresh-history
controls use `runtime-frozen-history/RANSPreservationProbe`; the primary frozen
solver and old replay executable remain untouched. The compact collector can
be reproduced against the retained immutable experiment with a new output path.

## Completion boundary and next work

This package qualifies diagnostic/admissibility corrections, independent
analytical weak identities, spatial approximation consistency and short
first-order behavior on the explicitly tested configurations. It supplies
bounded integration evidence, not exact discrete preservation of an arbitrary
polynomial reference. The unresolved unperturbed epsilon signals remain listed
as unresolved, even when the corresponding perturbed response has measured
first-order behavior.

A longer **staged neutral-persistence pilot** is justified only on a selected
h/p-resolved mesh with the measured adequate quadrature and monitored domain,
covariance, interface, trace-lag and drift budgets. Its purpose would be to
measure accumulation of the documented spatial and splitting errors. Before
calling that preservation evidence, advance the selected refined reference in
a short coarse/fine pair and confirm its evolved-state spatial residual/drift
trend. Do not extend the coarse arbitrary polynomial reference merely to hide
its initial residual. Long production campaigns, physical benchmarks, complete
temporal-order campaigns, terrain, true 3D, MPI and nonlinear/action solvers
remain outside this package; their unsupported combinations remain guarded.
