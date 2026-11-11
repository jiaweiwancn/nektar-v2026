# Native ownership correction and requalification

Generating-geometry ownership is corrected and verified through fresh native
integration. The six Natural short configurations and A95 Accepted neutral
qualification through T=.5 pass. All three bounded A95 responses now pass native
weak balance through .016. A93's common-w reduction and physical-budget accuracy
remain failed/unresolved independently, so a later response extension is not
eligible. No continuum equation or acceptance threshold is changed.

Starting revision is `5a1bc150b72df96bafb2f3cc8b75479064fd4bff`, with a clean
implementation worktree and matching remote. The unrelated main Mac checkout,
Linux primary checkout, CLS/WJ campaigns and every historical evidence directory
remain unchanged. Specifications retain hashes
`056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64`
and `17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40`.

The authoritative corrected-runtime run root is
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3`.
Its source/build/dependency symlinks name this package's isolated
`a93-a95-e94e95-ownership-20261111-export-v2` source/build/dependency area.
Frozen runtime, tmp and trajectories have fresh destinations. Earlier attempts
under `a93-a95-e94e95-ownership-20261111` and its `-export-v2` sibling retain
their frozen binaries, inputs, checkpoints and execution failures.
Numerical libraries use one thread, native execution is serial flat 2D, and
the owner-aware resource guard reserves startup headroom within 72 workers.

## Contract before changes and evolution

The legacy operator is labelled `LegacyCoordinateToleranceV1` in this report;
its checkpoints did not carry a separate ownership tag. The corrected numerical
operator is `GeneratedVolumeOwnerV1`, separate from the physical resolved
model descriptor. Its checkpoint and replay identities must reject old evolved
histories. Explicit read-only legacy-state evaluation records both identities
and cannot advance or resume a history.

Retain the physical equations, coefficients, kappa policies, finite inner and
zero high-Re outer molecular diffusion, full signed production and stress,
pressure convention, uncapped wall sources and one-sweep ordering. Eta and all
corrections remain zero. No clipping, snapping, panel removal, source cancellation,
reference resets or post-k epsilon repair. Keep the 3000-unknown and unsupported
route guards.

Retain each wall quadrature panel's endpoints before rounding its Gauss points.
Map that panel to its native column geometry, then use that owner consistently
in velocity/k and pressure rows. Retain outer quadrature's existing owner for
pressure and epsilon as well. Strict physical containment precedes trace-side
selection; a nearby interior point must not become a trace. Reject an explicit
wrong owner. Any local-coordinate endpoint arithmetic correction is restricted
to a physically contained point in a known owner and exposed diagnostically.
Coordinates, weights, cutoffs, physical breakpoints and source evaluations stay
unchanged. Cache keys retain owner, field, side and coordinate identities.

Establish production-path geometry regressions before evolution: committed
4/5 witness, permuted IDs, cutoff below/at/above a node, tiny representable panels,
exact-H sides and nearby interior points, all participating field owners,
partition/constant/volume/Stokes/free-row constant-flux/gauge checks. Repeat both
frozen Accepted/EP witnesses without altering their state, coordinates or weights.
Use unchanged analytical weak criteria. An actual perturbation need not have a
zero equilibrium residual.

Fresh corrected short cohorts cover six Natural model/wall configurations,
B5/B11 and dt=.0004/.0002, T=.0016/.016. Repeat selected Prescribed .41 controls
(A95 Accepted and A93 EP). Full native and fixed common dominant u/inner-k,
and resolved nonlinear-model w, retain .85 ratio, 1e-8 floor and 1e-10 absolute
uncertainty with the existing five-uncertainty resolution test. Require the
separate inner-only common-k check. Other modes remain visible. Neutral scaled
RMS limits are U/k/outer epsilon/Pi <=1e-3, w<=1e-4, mean wall/full interface
traction <=1e-3. Fieldwise temporal uncertainty remains
max(delta_t_B5,delta_t_B11)<=max(1e-7,.1*delta_h).

After each model's short gates, independently admit fresh A93/A95 Accepted
.002/.001 duration cohorts with .016/.1/.5 milestones and all four companions.
Only a physical-admissibility-preserving temporal failure can admit the already
permitted fresh .001/.0005 trial. Old evolution receipts do not qualify this
operator. A93 perturbation is outside this package, regardless of neutral outcome.

After complete corrected A95 Accepted .5 admission, repeat Accepted amplitudes
2e-5/1e-5 and A95 EP 2e-5 through .016 only, with independent corrected neutral
EP companions. Bind actual initialized perturbations and native histories
separately from neutral reference/boundaries. Compare all 17 baseline-subtracted
observables, with one gauge, spatial limit max(5e-7,.1*signal), temporal limit
max(1e-7,.05*signal), and half-amplitude limit max(5e-7,.1*half-signal).
Below-floor observables remain unresolved. No response extension is authorized.

Physical-budget reconstruction and accuracy are separate required scopes.
Repeat Accepted full-response volume48/64 and face16/32 controls on both spaces
and both dt values. Retain dependent epsilon-boundary storage, frozen/current
splitting and separately evaluated physical face fluxes, including natural-wall
gradient flux and cancellation of the internal interface in whole-patch totals.
Keep the existing .85 finite-defect ratio, max(5e-9,5*quadrature) floor and
max(1e-9,.1*spatial difference) timestep limits. Recompute primitives and bind
configuration/amplitude/horizon/members/native hashes and quadrature identities.
Reconstruction alone cannot qualify accuracy or admit a later response extension.

If A93 or physical-budget failures persist, diagnose signed assembled loads or
storage/source/physical faces before proposing at most one new spatial control
per independent mechanism, within the unchanged guards. Record its hypothesis
before execution. No open-ended mesh search or backend expansion.

## Historical status, retained unchanged

Old A95 Natural Accepted neutral .5 passes. Old A93 .5 common-w reduction fails.
Old A95 full/half Accepted and full EP .016 field comparisons pass while native
weak balance fails. Old physical-budget reconstruction passes and accuracy is
failed/unresolved. Four accepted-state/history failure witnesses and all 303
historical artifact hashes are checked unchanged. Eight historical epsilon-order
observations remain unresolved under their original criteria.

## Results

The production correction retains the panel endpoints in
`AWFProfileQuadraturePoint`. Each wall ray has an immutable native column/interval
catalog. `Refresh` chooses the generating panel's owner before evaluating its
rounded Gauss coordinates and maps that geometry into velocity/k/pressure.
Outer volume integration retains its existing generating geometry for all
participating fields, including epsilon. The evaluation cache retains field,
geometry, native element, exact side and coordinate bits. No source, cutoff,
point coordinate or quadrature-weight arithmetic changed.

`LocateElement` uses strict physical containment; only exact shared endpoints
can have multiple containing elements. Radial-side selection is separate from
interior containment. `ElementPoint` rejects a point outside the requested
native owner. A physically contained point's local-coordinate arithmetic may
exceed an endpoint by at most eight machine epsilons; any such local correction
is recorded with the raw local coordinate and an axis mask. No physical point
is snapped or accepted in a tolerance-expanded neighbour.

The strict rule exposed two additional coordinate searches in field export and
the driver's represented-reference error calculation. They now restrict native
coefficients through matching geometry, basis and signed CG maps, then use the
native backward transform. This is exact restriction of the existing polynomial,
not a change to the solution, its reference or a governing equation. Driver
errors still refer to the represented reference; the qualification's analytical
errors and independent weak oracle remain distinct.

The final native build succeeds at source
`7cb6757e06e6b1a619b01a98a2475dddb40a293d`, GCC 11.4 Release,
`-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS`.
Solver SHA256 is
`99828f96dcc95e0d529f315a136e72cec63b62f8e500e6c9017aa71712d4f914`;
preservation probe is
`56503e8e57d3b27b7583c05ef2a859a223d1c9be7ebb829bdac673f13be8d14b`.
Subsequent controller/report commits are recorded separately from compilation.

All focused geometry controls pass. The old committed interior witness selects
geometry 4; the corrected path selects its containing geometry 5. Permuted IDs,
cutoff below/at/above a node, exact-H/nearby samples and shared field ownership
pass. The smallest retained panel is `1.5612511283791264e-17`; no panel was
removed. Maximum constant-gradient, Stokes and gauge-volume defects across the
recorded controls are approximately `1.14e-14`, `5.87e-15` and `3.02e-14`,
against the predeclared `2e-12` primitive limits. Both frozen nonlinear Accepted
and EP witnesses pass the unchanged analytical weak tests on actual production
rows, with bitwise unchanged coordinates and weights. Their perturbed-state
spatial residuals are not required to vanish.

Same-version native restart, cache-on/off state and complete history equality,
and complete rollback pass. The real old-operator evolved checkpoint is rejected
with `checkpoint signature mismatch`. Ordinary CLS restart/cache/rollback
controls pass. All 11 targeted CTests pass, including native admissibility,
frozen boundary/interface controls, independent constitutive checks and the
required-scope negative controls.

Earlier output-stage failures are execution failures, retained separately from
physical or accuracy failures. A short-controller missing comparison-directory
failure is supplemented through explicitly identity-checked sealed reuse, not
by overwriting its failed job. A first debugger launch created 129 threads; the
owner-aware guard suspended it, it was terminated, and the reproduction used
one debugger worker. That debugger attempt is not claimed compliant with the
72-worker rule and supplies no qualification result. Scientific builds/runs
retain their independent worker censuses, CPU costs and suspension records.

The independent constitutive/wall kernel suite passes 1346 checks across all
three models, two walls and both coefficient policies. The selected fresh
Prescribed-kappa controls, A95 Accepted and A93 EP on B11 at dt=.0002 through
T=.0016, pass the model-specific common-test, field, domain and independent
stress/production checks. Their first audit correctly rejected a nonnested
default common space; the supplemented audit uses the explicitly recorded
model/policy-specific active-break common space and the same sealed trajectory.
This is selected short coverage, not a four-member Prescribed-duration campaign.

The A93 legacy T=.5 coefficients have also been evaluated read-only by the
corrected operator. State, quadrature coordinates and weights remain identical,
and no history is exported or advanced. Common-w B11/B5 changes only from
`1.01346325664` to `1.01346037467`, above the unchanged `.85` limit. The corrected
norms on those old coefficients are `1.08000740228e-7` and `1.09454470657e-7`.
This comparison does not qualify fresh evolution and shows why the remaining
A93 question must be tested independently from the geometry fix.

## Fresh corrected-operator neutral results

All six Natural configurations pass the B5/B11, dt=.0004/.0002 short matrix
through both T=.0016 and .016. All twelve complete comparison/audit bodies are
present. The supplemented export verifies each configuration, required time,
operator, source receipt hash and body identity. A pass-labelled control with
E94E95/EP's .016 observation removed fails with `missing_required_time:0.016`.

A95 Natural Accepted passes the complete .002/.001 four-member duration cohort
at .016, .1 and .5 (final theta=.1). Both final primary spatial audits, all four
complete scalar observers, both-space frozen quadrature and independent
constitutive checks pass. Each primary observer covers 251 unique accepted
states and 253 records; each half-step observer covers 501 unique states and
503 records. The two repeated restart boundaries have matching identities.

| A95 T=.5 scaled mass-dual residual | B5 | B11 | Ratio |
|---|---:|---:|---:|
| Native momentum-x | 1.126894e-6 | 3.764294e-7 | .334042 |
| Native inner k | 1.350592e-6 | 4.143729e-7 | .306808 |
| Native momentum-z | 7.182218e-7 | 2.693790e-7 | .375064 |
| Fixed-common momentum-x | 6.783273e-7 | 3.242799e-7 | .478058 |
| Fixed-common global k | 8.974764e-7 | 3.750535e-7 | .417898 |
| Fixed-common momentum-z | 2.501352e-7 | 1.158366e-7 | .463096 |
| Fixed-common inner-only k | 8.974749e-7 | 3.750529e-7 | .417898 |

These resolved ratios satisfy the unchanged .85 limit. Remaining equations
and new fine-space modes remain in the full receipts; no rate is assigned to
a below-floor quantity. Independent R/Q/P comparisons at the final primary
states have maximum absolute discrepancies of approximately 2.22e-16,
1.11e-16 and 4.17e-17. Full stresses include Ryy and production uses the raw
gradient, including isotropic divergence work.

| A95 T=.5 scaled reference weighted RMS | B5 | B11 |
|---|---:|---:|
| U | 1.947736e-5 | 7.642459e-7 |
| vertical velocity | 3.115986e-7 | 2.742108e-8 |
| k | 3.631649e-6 | 1.599723e-7 |
| outer total epsilon | 6.501799e-9 | 2.071413e-9 |
| augmented Pi, common gauge | 2.583396e-6 | 4.104976e-7 |

A95 full/half Accepted responses and the selected EP control are limited to
T=.016. EP requires its own corrected short neutral cohort before advancement.
A93's independent fresh duration cohort and the bounded response audits are
complete, with the separate outcomes below.

## Bounded full-response result and independent budgets

A95 Accepted amplitude 2e-5 passes all four field/response comparisons at T=.016,
and its complete native audit now passes on both spatial members. This includes
the analytical weak tests that failed under the legacy operator, all four scalar
observers, frozen quadrature and independent constitutive checks. Fourteen of
the seventeen response observables are resolved; k, outer epsilon and Qyy are
below their declared signal floors. Planar A95 Qyy is exactly zero here, while
Ryy remains explicitly retained and has a resolved response.

For example, primary Qxz response RMS is approximately 1.427917e-4; its spatial
difference is 3.998991e-8 and maximum timestep difference 7.045941e-7, within
the 7.151276e-6 temporal limit. Primary augmented-pressure response RMS is
3.773466e-4. These are baseline-subtracted responses on fixed physical scales,
not temporal-order estimates or interface-traction-jump errors. Neutral Qxz=0
does not mean the full neutral outer Rxz vanishes.

All sixteen required full-response physical-budget reconstructions pass on the
corrected runtime, including volume orders 48/64 and face orders 16/32. Their
finite-defect accuracy still fails the original derived classifier. The exact
binding checks configuration, amplitude, horizon, four native requests/reports,
operator and integration controls, and recomputes all primitive classifications.
Reconstruction alone cannot qualify accuracy or admit a response extension.

| Frozen-stage signed scaled defect | B5 | B11 | Spatial / timestep status |
|---|---:|---:|---|
| Perturbed inner k | 1.227549e-5 | 1.232417e-5 | failed / unresolved |
| Perturbed outer k | 1.426296e-5 | 1.426437e-5 | failed / unresolved |
| Perturbed outer epsilon | 1.389703e-5 | 1.389741e-5 | failed / unresolved |
| Perturbed whole-patch k | 1.323038e-5 | 1.325075e-5 | failed / unresolved |
| Neutral inner k | -1.702987e-6 | -1.654268e-6 | failed / unresolved |
| Neutral outer k | -3.102471e-9 | -1.680538e-9 | below floor / passed |
| Neutral outer epsilon | 1.464027e-7 | 1.468634e-7 | failed / unresolved |
| Neutral whole-patch k | -8.224477e-7 | -8.020719e-7 | failed / passed |

The near equality with the archived defects is an observed result, not reuse of
the old evolution. All current states and histories are fresh. Physical budgets
retain actual consecutive-state storage, dependent epsilon-boundary storage,
frozen/current source and coefficient splitting, separately integrated physical
advective/diffusive face fluxes and nonzero reconstructed natural-wall gradient
flux. Inner/outer interface terms remain visible and are excluded from the
whole-patch external total. No prescribed-row reaction replaces a physical flux.
The epsilon constant-test identity on this one-column prescribed-side patch is
labelled a reaction identity, not an independent conservation theorem.

## Remaining A93 spatial gate

All four fresh A93 Accepted trajectories reach T=.5 with passing physical-domain,
field, traction, source/storage, coverage, constitutive and quadrature checks.
The complete duration qualification fails only fixed-common vertical momentum.
This is a post-step qualification failure in accepted physical states, not a
rejected timestep. No further A93 evolution or perturbation follows it.

| A93 T=.5 reduction | B5 | B11 | Ratio / status |
|---|---:|---:|---|
| Native momentum-x | 6.090692e-7 | 3.847890e-7 | .631766 / passed |
| Native inner k | 9.697020e-7 | 4.744441e-7 | .489268 / passed |
| Native momentum-z | 5.232139e-7 | 2.866422e-7 | .547849 / passed |
| Fixed-common momentum-x | 3.970654e-7 | 3.305460e-7 | .832472 / passed |
| Fixed-common global k | 7.111434e-7 | 4.338758e-7 | .610110 / passed |
| Fixed-common momentum-z | 1.080002e-7 | 1.094542e-7 | **1.013463 / failed** |
| Fixed-common inner-only k | 7.111421e-7 | 4.338754e-7 | .610111 / passed |

The half-timestep common-w norms are 1.102536e-7 and 1.113365e-7; their ratio
also fails .85. These values are above the unchanged 1e-8 residual floor.
The signed spatial load difference is 1.734775e-7. Its timestep-load differences
are 1.377019e-8 and 7.925570e-9, below the separately declared 1.734775e-8
uncertainty limit. Thus the retained failure is resolved by that diagnostic;
halving the timestep does not remove it.

Actual assembled pressure and constitutive-stress pieces each have common-dual
norm about .692438 and largely cancel. Their individual norms are not numerical
defects. The globally constant total-flux free norms are only 1.66e-13 and
2.28e-13, far below the remaining load. Signed pressure/stress/advection projections,
full term reconstruction, constant/varying-flux element buckets and actual row
supports are retained for all four states. Dominant common rows span
z=.04.. .5, particularly across roughness height .1. This supports a localized
spatial/coupling investigation; it does not identify the geometric ownership
defect as the cause of this remaining failure.

## Bounded budget direction controls

Before either control was initialized, two separate hypotheses were recorded:
normal endpoint-cell resolution (split only 0..0.00125 and 3.5..4), and streamwise
degree at fixed vertical space. Their physical-integral comparisons retain the
original floors and uncertainty limits and cannot overwrite the required B5/B11
failure. No field arrays on different meshes are zipped or truncated.

The normal-endpoint control has 2589 actual native unknowns, unchanged p5/p3
spaces and all original active breaks. Four independently initialized histories
(neutral/2e-5 perturbation, dt=.002/.001) reach .016 with passing native accuracy
and domain checks. Its eight frozen physical-budget controls reconstruct and
resolve quadrature. Neutral inner-k, outer-epsilon and whole-k defect magnitudes
fall to .096707, .065620 and .064305 of B5 respectively; their timestep checks
pass. Neutral outer k remains below the original floor.

Perturbed inner/outer-k, epsilon and whole-k defect ratios are respectively
1.12531, .993495, .964744 and 1.05463; none qualifies. Their timestep accuracy
also remains unresolved. The baseline-subtracted inner-k defect changes only
from 1.39784785e-5 to 1.39784287e-5, about 5e-11. Removing a cancelling baseline
endpoint error therefore does not remove the distinct response-related defect.
The signed storage/source/region/physical-face terms are retained, rather than
attributing every defect to one face or to source splitting.

The streamwise p6/z-p5 control predicts 2769 unknowns with pressure p4/p3, but
the native initialization rejects it with
`Global wall field: requires equal-order Modified_A tensor bases`.
The actual count is unverified because this guard fires before patch assembly.
There is no evolved state or history. Its complete input, native failure packet
and log are archived. This is an unsupported spatial-basis capability, not a
constitutive-domain failure or a successful expected negative control. The
original failed execution remains unchanged. No guard is relaxed and no
replacement mesh search follows. Same-vertical-space two-column B5/B11 would
require 4295/5466 unknowns, exceeding 3000. A clean isolated streamwise-resolution
experiment remains blocked within the present supported spaces.

## Response coverage and interpretation

| Corrected-operator scope | Maximum T | Outcome |
|---|---:|---|
| Six Natural model/wall short cohorts, B5/B11, .0004/.0002 | .016 | All twelve required time/configuration bodies pass |
| A95 Accepted neutral, B5/B11, .002/.001 | .5 | Complete duration qualification passes |
| A93 Accepted neutral, B5/B11, .002/.001 | .5 | Native integration completes; common-w reduction fails |
| A95 EP neutral companions, B5/B11, .002/.001 | .016 | Pass |
| A95 Accepted responses, amplitudes 2e-5 and 1e-5 | .016 | Field uncertainty and complete native audits pass |
| A95 EP response, amplitude 2e-5 | .016 | Field uncertainty and complete native audits pass |
| Required Accepted full-response physical budgets | .016 | Reconstruction passes; accuracy fails/unresolved |
| Selected Prescribed .41 controls: A95 Accepted, A93 EP | .0016 | B11/.0002 checks pass; no duration cohort |

B5 and B11 retain 2355 and 2997 unknowns respectively, p5 velocity/scalars with
compatible p3 pressure. Their exact interval lists, active breaks, quadrature,
space hashes and configuration-specific common tests are bound in each request.
The normal-endpoint diagnostic is a separate 2589-unknown space, not a replacement
qualification pair.

All three A95 Natural cohorts pass through .016: Accepted amplitudes 2e-5 and
1e-5, and EP amplitude 2e-5. Each has both spaces, both timesteps, actual perturbed
initial-state/zero-step-history binding, complete native weak/spatial/scalar and
frozen-quadrature audits, independent full-tensor/source comparisons and zero
correction activity. EP has its own fresh .002/.001 neutral companions; no old
or Accepted-wall history admits it. All 68 half-amplitude checks pass. The common
physical-grid EP/Accepted comparison retains all 17 observables, explicit
interface sides and one pressure gauge; each wall's uncertainty checks pass.
Different wall-law responses are not required to be equal.

These results close the previous native analytical weak-balance response failures.
They do not close physical-budget accuracy, establish temporal order, qualify
A93 perturbations, extend EP neutral duration to .5, or qualify any response
beyond .016. The eight historical unresolved epsilon-order observations remain
unchanged. New timestep pairs bound uncertainty only.

## Provenance, evidence and cost

`NativeOwnershipRequiredScopeV1` keeps ownership, history/cache/CLS controls,
short-matrix coverage, each neutral model, response fields, native weak balance,
physical reconstruction and physical accuracy separate. Its overall status is
failed/unresolved because the independent A93 and budget gates do not pass.
`response_extension_eligible` is false. Ten explicit required-scope controls
pass, including missing/foreign/stale accuracy, wrong horizon, missing companions,
pass-labelled primitive mismatches, old-operator inputs and a command-line
reconstruction-pass/accuracy-fail fixture that exits 1. Actual old native history
rejection is separately exercised, not inferred from the fixture.

Exact completed-request verification returns 0 for the complete short matrix
and A95 final audit, and 1 for A93's failed final audit. Receipt hashes remain
unchanged before/after these checks. Full numerical comparison bodies and
classifications are exported, including required times and consumed initial,
reference, boundary and history identities. All-state/stage extrema and cutoff
motion retain their native stage/time/side/geometry metadata; supplementary
fixed-grid scalar observers retain both minima and maxima with coordinates.
Final-state diagnostics remain separate. Native finite-inventory eligibility is
not replaced by the supplementary sampling grid.

The final source audit matches all 209 recorded compiled C++/header hashes to
the frozen 7cb6757 snapshot. Later changes are isolated qualification/reporting
drivers. Requests record exact scientific script hashes independently of Git
revision. Frozen executable identities are:

| Artifact | SHA256 |
|---|---|
| IncNavierStokesSolver | `99828f96dcc95e0d529f315a136e72cec63b62f8e500e6c9017aa71712d4f914` |
| RANSPreservationProbe | `56503e8e57d3b27b7583c05ef2a859a223d1c9be7ebb829bdac673f13be8d14b` |
| RANSNewModelsProbe | `fd7d3346e68fd67138f4755b3f5f6c0513927a2545651d63213d98d1030d85b0` |
| RANSNativeDiagnosticChecks | `4d62a4b37b9a36c326bec5826033eacad53d98b873c54a284e5f6e4dd12aca20` |

Actual `/proc` loader maps include native solvers and probes, plugins and thread
environments. The legacy and corrected probe executables carry their distinct
ownership implementations statically; all 34 compared copied shared dependencies
are byte-identical to the legacy private copies. No shared solver installation
or campaign library is replaced. Runtime aliases ending in v3/v5 are recorded
private copies of the corrected snapshot, not use of the historical operator.
Linux Python 3.10.12 uses NumPy 1.21.5 and SciPy 1.8.0; compiler/build flags and
SDK/library hashes accompany the source manifests. MPI-enabled dependencies are
used in rank-one execution with numerical-library threads explicitly one.

Ninety-four completed numerical segments contain 4696 accepted steps, including
the bounded endpoint diagnostic control. Their instrumented exclusive process
CPU totals 31812.1444 seconds (8.84 CPU hours). Diagnostic probes, compilation,
orchestration and retired attempts are additional and retain separate receipts.
`domain.compute` alone accounts for 15775.861 seconds; flow assembly 3901.355,
residual work 2528.732, complete audits 2503.927 and physical audits 1822.626.
Dense factorization plus triangular solves total 456.552 seconds, about 1.44%.
Output/checkpoint scopes and per-job user/system CPU remain in the machine data.
No backend or state-dependent validity-cache optimization is introduced.

Guard-recorded solver suspension sums to approximately 1924.969 worker-seconds,
with no unclosed solver suspension. This is not campaign elapsed time; parallel
intervals can overlap. Gross wall time includes waits and suspensions. Scientific
guard census peaks at 62 active workers, with startup headroom inside 72. The
separate earlier debugger exception remains explicitly retained. The final owned
numerical census is empty; only this experiment's guard is stopped with SIGTERM,
avoiding any signal to stale paused IDs or external campaigns.

The adjacent `a93-a95-e94e95-awf-native-ownership-requalification-evidence` contains
`decision.json`, complete compressed result bodies, all-state/stage receipts,
source/runtime/input/history identities, failure replay requests and native
checkpoint copies, geometry/unsupported-guard fixtures, and exact Linux controller
sources. `artifact-index.json` maps heavy retained packets/reports to their Linux
paths and hashes. `SHA256SUMS` covers the compact export. All 303 historical
artifact hashes and four archived failure witnesses remain unchanged; no old
success is promoted to qualification of this new numerical operator.

The final collector runs at `1333017f1dd4cb7fd4a34393ab42d44d5487d0cf` and
exports 203 complete result bodies and 94 all-state/stage segment receipts.
Independent export validation checks all payload hashes, all twelve short
bodies and the 209 unchanged compiled-source hashes. The actual aggregate CLI
returns 1 with exactly `A93_neutral` and `physical_budget_accuracy` unmet;
`collection-execution.json` and `export-validation.json` preserve that result.
Its nonzero exit is the scientific outcome, not an expected-negative-control
success. Exporting takes 76.071 wall seconds and 56.365 child CPU seconds,
separate from numerical evolution.

## Reproduction

Use the preserved Linux root and private runtime, not an active CLS installation:

```bash
R='/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3'
S="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 TMPDIR="$R/tmp"
export LD_LIBRARY_PATH="$R/runtime/lib"
```

These exact-request checks verify completed evidence without advancing a state:

```bash
python3 -B "$S/run_rans_ownership_short.py" --root "$R" \
  --geometry-receipt "$R/geometric-regressions-v1/results.json" \
  --concurrency 2 --output-label ownership-short-qualification-v2 \
  --audit-prefix short-audits-v2 \
  --reuse-executions "$R/ownership-short-qualification" --resume
python3 -B "$S/finish_rans_balanced_duration.py" --root "$R" \
  --admission "$R/admissions/A95-N-Accepted.json" --dt .002 --time .5 --resume
# Expected exit 1: the original common-w gate remains failed.
python3 -B "$S/finish_rans_balanced_duration.py" --root "$R" \
  --admission "$R/admissions/A93-N-Accepted.json" --dt .002 --time .5 --resume
python3 -B "$S/diagnose_rans_duration_residual.py" --root "$R" \
  --admission "$R/admissions/A93-N-Accepted.json" \
  --cohort A93-N-Accepted-duration-0.002 --time .5 \
  --output "$R/duration-residual-diagnosis/A93-GeneratedOwnerV1-0.5" --resume
```

The last diagnostic also returns 1 while preserving its completed decomposition.
Recompute budget accuracy in a new destination (expected exit 1):

```bash
python3 -B "$S/classify_rans_physical_budget_accuracy.py" \
  --source "$R/physical-budget-analyses/A95-N-Accepted-response-2e-05-dt-0.002-0.016/results.json" \
  --output "$R/reproduction-budget-accuracy-new"
```

For native reproduction, the complete argv, immutable requests, configurations,
zero-step bindings and continuation controllers are preserved in the export and
Linux run root. New numerical runs require a new destination and a fresh worker
census under the recorded owner-aware guard; only exact validated resume may
reuse completed artifacts. Compressed failure requests can be expanded beside
their copied XML and passed directly to the frozen preservation probe for
read-only replay. They retain the trajectory's consumed boundary map. Native
checkpoints are copied byte-for-byte, never re-signed or repurposed across
operator, mesh, timestep or method identities.

## Next bounded gate

The ownership repair is justified by independent geometry identities, unchanged
frozen point/weight witnesses and fresh integration. It does not resolve A93's
remaining common-w convergence or the response-related physical-budget defect.
No longer response is justified by this package.

For A93, a bounded next experiment can redistribute one nonessential refinement
within the 2997-unknown p5/p3 fine mesh toward the measured roughness-linked row
supports, using the retained varying-flux buckets to select the donor interval.
It must preserve active breaks, use fresh histories, and retain the same .85
common-w and inner-only k gates. No coefficient or source alteration is indicated.

For budgets, first qualify an affordable streamwise-resolution capability. The
2769-unknown anisotropic proposal is blocked by the current basis guard; a
two-column p5 counterpart exceeds 3000. Any basis/backend extension therefore
needs a separate implementation and geometry/restart qualification package,
followed by the same independent physical-budget accuracy test. The baseline
endpoint improvement alone cannot admit a response extension.

After those gates, a multi-column neutral inflow/outflow test must precede curved
terrain. Current success concerns a prescribed-reference bounded patch only.
Curved geometry, true 3D, MPI, general nonlinear/action solvers, open boundaries,
A93 perturbations, long-time stability, global realizability and physical
benchmark validation are not attempted or claimed here.
