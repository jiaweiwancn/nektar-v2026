# M8 compact proof and diagnostic evidence

## Material Passport

- Origin: approved M8 model proposal and academic-research-suite experiment-agent workflow.
- Date: 2026-11-03.
- Mathematical/implementation verification: **completed for the receipts listed here**.
- Broad robust-model qualification: **FAILED for the screened variants on the declared state box**.
- Physical terrain accuracy: **not qualified by these experiments**.
- Authoritative numerical environment: Linux through Tailscale, `jay@100.105.247.110`.

Passing a verification receipt means that its identities, numerical comparisons,
or explicit failure witnesses were reproduced. It does not turn a failed model
property into a passing qualification. The milestone report owns the complete
native/regression campaign and final integration status; this index describes
the compact mathematical, native, reduced-mode and state-coverage evidence.

## Files and interpretation

| Evidence | What it establishes |
|---|---|
| [Fresh kernel qualification](proof-kernel-fresh.json) and [runner](proof-kernel-fresh-runner.json) | 37,310 independent checks plus 27,907 native probe checks; explicit common-cubic scaling, consistent modified production, capped-source identities and one-sided actions, required invalid-input rejection. |
| [Kernel replay](proof-kernel-replay.json) and [native provenance](proof-native-replay-provenance.json) | 37,309 checks with the later interval-proof code and the same stored native report. One native-command exit check is absent because this is a replay. Original native execution `07363199` is distinguished from the AA compiled campaign and `11aa94c0` proof source. |
| [Transverse qualification](transverse-validation.json) and [frozen plan](transverse-plan.json) | 305,254 checks of the reduced Fourier/native-CG three-component transverse subsystem, with independent stress, weak matrices, native histories, energy, covariance and refinement oracles. Local decay is qualified; coupled/terrain stability is not. |
| [Joint local proof, AA](proof-local-joint-aa.json) | Original local momentum/covariance certificate used by the recorded native state-coverage audits. |
| [Joint local proof, 11aa](proof-local-joint-11aa.json) and [runner](proof-correlated-runner.json) | Sharper directed-rounding coefficient enclosures, the same declared local boxes, fixed-gamma covariance proofs, and an explicit outside-envelope control. |
| [Conservative wide bound](proof-conservative-wide.json) | Sufficient gamma values from conservative interval bounds. These values are **not necessary or minimal**; a failure of an overdiffusive sufficient choice alone would not rule out every gamma. |
| [Fixed-gamma feasibility proof](proof-fixed-gamma-feasibility.json) and [plan](proof-fixed-gamma-plan.json) | Separate exact rational witnesses impose contradictory **necessary** lower/upper gamma constraints. This stronger result rules out the specified fixed scalar gamma law on the declared broad box. |
| [Invariant-response analytic controls](proof-invariant-controls.json), [runner](proof-invariant-controls-runner.json), [independent screen](proof-invariant-screen.json) and [plan](proof-invariant-plan.json) | Explicit directed-interval constitutive/covariance failures of the two optional C(q) research variants, with all C′ and β′ terms included. No V3 production model was adopted. |
| [p3 reference coverage](coverage-reference-p3.json), [perturbed coverage](coverage-perturbation.json) and [runner](coverage-runner.json) | Independent reconstruction at the actual declared outer source, mass, reconstruction, interface and transfer evaluation sets. These p3 cases leave the local certificate box. |
| [p5 coverage](coverage-reference-p5.json), [p7 coverage](coverage-reference-p7.json), [refinement runner](reference-refinement-runner.json), [p5 energy](reference-p5-energy.json), [p7 energy](reference-p7-energy.json) | Actual zero-step reference projections with independent energy checks. p7 covers the finite outer certificate sets and has positive sampled inner covariance; this is an initial-state result. |
| [Runtime epoch summary](runtime-epoch-summary.json) and [build runner](runtime-build-runner.json) | Compact source/runtime binding for the producing build epoch. The complete fingerprint remains on Linux; this summary does not replace final source/runtime equivalence checks. |

## Reduced transverse experiment

The figure selects the actual `p7-n8-xi4`, IMEX2, amplitude `1e-4` rows at
`dt=0.00015625 s`, `T=0.05 s`. The disturbance is a reduced Fourier/native-CG
experiment supported above the wall layer, not a nonlinear terrain trajectory.

| Model in the reduced experiment | Final energy ratio |
|---|---:|
| Original CLS, reused M7 control | 8.62184899621 |
| Common cubic β=0.4, Γ=0 | 0.589289214077 |
| Common cubic β=0.45, Γ=0 | 0.734610537997 |
| Common cubic β=0.4, Γ=0.25 | 0.333442685761 |

Here Γ multiplies the prescribed neutral eddy viscosity in the persistent
reference-preserving stress increment. Γ=0.25 is a reduced-mode sensitivity
case. It is distinct from the Γ=1 native reference-projection cases and from
the fixed-gamma feasibility problem. The original control is the preserved
[M7 receipt](../m7/transverse-validation.json), not a new M8 solver execution.

## Local certificates and broad fixed-gamma failure

Use the upstream orthonormal frame and normalization
`A*=A/G0`, `k*=Qstar*k/K0`, `epsilon*=epsilon/epsilon0`, where
`G0=u_star/(kappa*r)`, `K0=Qstar*u_star^2` and
`epsilon0=u_star^3/(kappa*r)`. The tangent has physical units νt,0 and
the Reynolds covariance has physical scale `u_star^2`.

Both limiting max-invariant branches are covered by the interval certificates.
The following 11aa bounds use `k/K0, epsilon/epsilon0` in `[.99,1.01]` and
entrywise gradient distance from `e_x tensor e_z`:

| Fixed Γ | Gradient radius | β=0.4 lower principal bound / νt,0 | β=0.45 lower principal bound / νt,0 | Covariance |
|---|---:|---:|---:|---|
| 1 | .001 | .471272810698 | .359364129975 | Strictly positive throughout the declared box by interval LDL |
| 3 | .01 | .547381238192 | .252563404664 | Strictly positive throughout the declared box by interval LDL |

These are local constitutive properties with frozen k and epsilon. They do not
prove state-box invariance, outer turbulence positivity, a high-order discrete
maximum principle, long-time coupled stability, or terrain accuracy. Covariance
positivity for one fixed Γ does not automatically extend to larger Γ.

The wider feasibility box has entrywise gradient radius `.25`, exactly
trace-free witness gradients, and `k/K0, epsilon/epsilon0` in `[.5,1.5]`.
The selected exact rational witnesses imply these necessary constraints:

| β | Momentum requires Γ at least approximately | Covariance requires Γ at most approximately |
|---|---:|---:|
| .4 | 2.44640920758313 | .177525831861587 |
| .45 | 2.27794509320359 | .156117116286741 |

The intersections are empty. The actual directed interval endpoints, rational
states and directions, branch checks and gamma-shift identities are retained
in `proof-fixed-gamma-feasibility.json`. The plotted intervals are necessary
conditions, not asserted sufficient or safe regions. This result concerns
the specified fixed scalar stress increment on this box, not every possible
constitutive regularization or global model impossibility.

## Optional invariant-response alternatives

Both research alternatives retain the original C(q) up to `qc=3.3` and use
the C1 power tail `C(qc)*(q/qc)^(-alpha)` above it, with
`alpha=-qc*Cprime(qc)/C(qc)`. Candidate A keeps common β constant; B uses
`beta0*[qc*C(qc)/(q*Ctail(q))]^2` above the join. B is the explicitly
specified research completion of the adaptive-coefficient direction. Neither
alternative was adopted in production.

Candidate A has an explicit adverse transverse coefficient at q=8. Candidate B
repairs that selected direction, but the independent strict-active screen
finds other negative incompressible principal forms for both beta0 values.
These forms fail the declared constitutive coercivity gate; they are not
new native coupled PDE trajectories or nonlinear instability theorems.

There is also a simple common covariance witness for both alternatives:
`A02=5/4`, other gradient entries zero, `k=19/4`, `epsilon=1/2`,
`q=95/8`, and covariance direction `(1,0,5)/sqrt(26)`. It lies within the
declared broad box and gives variance `-0.2832467961122523`. The screen's
boundary control with `k=1.5*Qstar` is a different state and gives approximately
`-0.310052142763282`; the two results are consistent.

## Initial projection and finite coverage

These results use β=.4, Γ=1 and the same physical model inputs:

| Degree | Maximum normalized gradient-entry deviation | Minimum inner covariance at source quadrature | Finite outer joint-certificate coverage |
|---|---:|---:|---|
| p3 | .0642009880560146 | -.177666349051788 | Outside |
| p5 | .00264835174984812 | -.0198538557236063 | Outside |
| p7 | .000102433249649092 | +.0806939841621756 | Covered at the declared finite sets |

The p3 figure row selects t=0 from the stored anchored control. Its subsequent
control coefficients are exactly identical, so the stored aggregate coverage
extrema equal its initial projection. The p5/p7 sessions have zero steps.
This separates resolution of the projected reference from the broad
constitutive failures. A positive p7 initial sample does not qualify an evolved
trajectory or guarantee the minimum over every point of the polynomial.

## Figures, curation and reproduction

- [Three-panel PNG](m8-qualification.png) and [vector PDF](m8-qualification.pdf).
- [Exact plotted values and input hashes](figure-data.json).
- [Figure receipt](figure-receipt.json) records the plot-source hash, producing
  data commits, versions, one-worker environment, output hashes and actual command.
- [Plot/curation source](plot-m8.py) reads actual immutable JSONs. It executes no solver.
- [Curation index](curation-index.json) records each full original path, byte
  count and SHA256 alongside the curated-file hash and declared transformation.

The final scientific figures are generated on Linux under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/figures-002`.
The initial layout attempt remains under `figures-001`. The final compact
archive is `curated-proof-002`. All rendering and mathematical proof jobs use
one worker with numerical library threads one; the figure run is CPU 48.

To render the exact extracted data again on Linux, use a fresh output path:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  taskset -c 48 python3 -B plot-m8.py --data figure-data.json --output FRESH_FIGURE_DIRECTORY
```

Large native fields, complete histories, solver logs and full 9x9 interval
tangent arrays remain in the Linux run directories. The three local/wide
proof JSONs remove only the recursively named `tangent_intervals` arrays and
add explicit curation metadata; their exact scalar bounds, null-form weights,
interval LDL pivots, states and covariance witnesses are retained. Other
listed proof/coverage receipts are copied byte for byte. The runtime summary
declares its omitted full manifests and the authoritative original hash.

## Native curved operator verification

The [curved curation index](curved-curation-index.json) records exact Linux
paths, bytes and SHA256 hashes for the byte-for-byte compact copies below.
Raw probe coefficients, meshes, directional vectors and solver logs remain
on Linux under `m8/verification/`. These are **nine distinct selected
probe/oracle cases**: three Original controls, one CubicScaled case and five
Regularized cases. They are not nine developed CFD trajectories.

| Actual case ID | Selected model and scope | Independent operator audit |
|---|---|---|
| curved-001/original-variable | Original, 2D Gaussian and variable H | [Passed](curved-001-original-variable.json) |
| curved-001/scaled-variable | β=.4, 2D Gaussian and variable H | [Passed](curved-001-scaled-variable.json) |
| curved-001/regularized-variable | β=.4, Γ=1, 2D Gaussian and variable H | [Passed](curved-001-regularized-variable.json) |
| curved-001/original-transition | Original, 2D cosine feet and broad 100/200 sensors | [Weak/action comparison passed](curved-001-original-transition-broad.json) |
| curved-001/regularized-transition | β=.4, Γ=1, same broad-sensor feet fixture | [Passed](curved-001-regularized-transition-broad.json) |
| curved-001/regularized-surface3d | β=.4, Γ=1, 3D Gaussian and variable H | [Passed](curved-001-regularized-surface3d.json) |
| curved-002/original-transition-active | Original, 2D feet and active 1/30 sensor fixture | [Passed](curved-002-original-transition-active.json) |
| curved-002/regularized-transition-active | β=.4, Γ=1, active 1/30 sensor fixture | [Passed](curved-002-regularized-transition-active.json) |
| curved-003/regularized-plane-rotated | β=.4, Γ=1, incoming Plane slope .25 | [Passed](curved-003-regularized-plane-rotated.json) |

These cases use probe source `6284787bf4c436203f4404498ba3dfb68aa144ef`.
The isolated Release MPI build compiles only `RANSCurvedWallProbe`, using
16 compiler workers. Numerical probe and oracle commands run sequentially
on CPU 50 with one numerical library thread. The primary AA solver and
geometry overlay remain unchanged. [Loader bindings](curved-probe-ldd.txt)
and all before/after preservation checks are retained in the
[001](curved-001-runner.json), [002](curved-002-runner.json) and
[003](curved-003-runner.json) runner receipts. The exact executed recipes
are [001](run-curved-001.py), [002](run-curved-002.py) and
[003](run-curved-003.py); they refuse to overwrite their historical run
directories. Reproduction needs fresh output paths and a new runner hash.

All nine selected weak/operator audits pass 212 recorded assertions. The
largest weak-load component discrepancy is `6.67e-15`; complete 2D coupled
velocity/k/epsilon/pressure directional comparisons are below `6.0e-11`.
The selected constitutive increment uses the same stress in momentum and
turbulence production. The independent source oracle explicitly verifies
that every evaluated ray remains in the open inactive repair plateau before
reusing the earlier source law; this includes 12,368 ray points for a 2D
upstream state and 20,480 for the 3D state, with minimum k/K above .95.

The [Original M5 control](curved-original-m5-control.json) passes 502 checks.
The [active-sensor Original M6 control](curved-original-m6-control-active.json)
passes 50 weak and 320 additional checks, including nonzero omitted-chi
derivative controls. The first broad 100/200 fixture did not activate that
legacy sensitivity control. Its [failure log](curved-001-original-transition-broad-legacy.log)
is preserved; fresh 1/30 fixtures use the existing M6 operator-test recipe
and unchanged acceptance tolerances. This was a fixture-sensitivity issue,
separate from the successful weak-load and action comparisons in that attempt.

The [Plane frame audit](curved-003-reference-frame.json) independently
reconstructs the incoming frame from the source XML and analytic geometry.
The flow is `(0.9701425,0,0.2425356)` and the normal is
`(-0.2425356,0,0.9701425)`; both metadata comparisons have zero error.

The [integration/remaining-gates record](integration-remaining-gates.json)
distinguishes these implementation checks from broad robust qualification.
The probe fields are bounded **unsolved generic states**. Their all-mode
operator audits establish neither a genuine curved-gamma manufactured
solution nor long-time or physical terrain qualification. The genuine
curved-gamma MMS gate remains unestablished.

## Corrected flat-gamma initialization and checkpoint compatibility

The corrected source `ed34052fe4dd967e963c22272a0ec0888cf35b70` passes
[seven fresh flat-gamma MMS cases](mms-initialization-qualification.json):
flat, p4, h2, dt2, MPI 2, checkpoint and restart. Every solver case executes
afresh with Regularized beta .4 and GammaScale 1. The driver records 270
checks; the largest independently integrated analytic field error is
`3.71e-16`. MPI and restart full-field differences are exactly zero.
Complete native-history mass-norm maxima are `1.57e-15` for MPI and zero
for the fresh restart. The checkpoint state also passes the independent
analytic MMS verification.

The [exact compatibility receipt](mms-checkpoint-compatibility.json) retains
the complete old/new signature strings, their SHA256 hashes, binary checkpoint
hashes/checksums, restart report hashes and canonical command hashes.
The stale AA flat Functions+Capped checkpoint is rejected before an accepted
restart. Its signature differs from the corrected signature only by
`:capped_initial_trace=1`. The unaffected old Neutral CubicScaled checkpoint
remains compatible and retains exactly the same signature. Its continuation
matches the old uninterrupted trajectory to `1.01e-14` in full-field norms
and `3.08e-17` over every native history slot, below the frozen `1e-8` limit.
This is an expected versioning rejection plus a positive compatibility control.

The [runner receipt](mms-initialization-runner.json) records source, solver,
commands, timeouts and preservation checks. At most two MPI ranks run with
one numerical library thread each; actual [MPI bindings](mms-initialization-mpi-bindings.log)
place ranks on socket 0 cores 0 and 1. Python is constrained to CPUs 40–43;
MPI may repin its ranks. Full before/after fingerprints retain 3,920 source
inputs and 100 runtime images, with unchanged source/images, loader bindings,
packages and tools. Binding reports intentionally differ from the build-002
environment. The separate read-only audit explicitly initializes the selected
oracle's VECLIB thread variable before its own pin and records that difference.
Actual duplicated/trailing-empty LD paths are preserved, and loader bindings
and all transitive library bytes are checked from all nine solver case directories.

The first [MPI launch attempt](mms-initialization-attempt-001.json) stopped
because its CPU-list environment conflicted with the driver's mapping flag;
its [launcher log](mms-initialization-launcher-failure.log) and four successful
serial prefixes remain archived. All seven cases were then executed afresh.
[Environment-audit attempts](mms-initialization-environment-audits.json)
also preserve the two postprocessing guards and the successful final audit.
The [small affinity diagnostic](mms-affinity-diagnostic.json) records its exact
argv and known environment overrides without claiming an inherited-environment
snapshot. The [curation index](mms-curation-index.json) pins all exact copies;
the executed recipes are [first runner](run-initialization-mms-001.py),
[final runner](run-initialization-mms-002.py) and
[final compatibility audit](audit-initialization-compatibility-003.py).
These initialization and flat-MMS results do not establish the separate
curved-gamma, broad admissibility, long-time or physical terrain gates.

## Native campaigns and source identity

The [native curation index](native-curation-index.json) binds each compact
file below to its original Linux bytes and exact declared transformation.
The [collector](curate-native.py) runs on Linux and refuses to archive an
unfinished campaign as complete. Native coefficient vectors, field files,
checkpoints and complete histories stay in the experiment directories.

| Receipt | Source-bound result |
|---|---|
| [CubicScaled native](native-scaled-validation.json) | 1,219 checks, 18 verified cases, IMEX2 difference ratio 4.00066029. Actual AA solver with corrected 628 audits. |
| [Regularized native](native-regularized-validation.json) | 1,214 checks, 17 verified cases, IMEX2 difference ratio 4.00826716. Same beta=.4, GammaScale=1 and AA solver. |
| [Native completion runner](native-completion-runner.json), [before](native-completion-runtime-before.json) and [after](native-completion-runtime-after.json) | Two concurrent jobs, each at most MPI2/one numerical thread per rank; all 100 runtime images and actual bindings unchanged. |
| [Original flat MMS completion](mms-completion-validation.json) and [protocol](mms-completion-preregistration.json) | 309 checks, seven cases, exact all-field MPI/restart agreement. This predates the startup correction; the fresh corrected-source matrix is separate. |
| [Unforced flow diagnostic](flow-validation.json), [runner](flow-runner.json), [before](flow-runtime-before.json), [after](flow-runtime-after.json) and [protocol](flow-preregistration.json) | Stopped at t=.055 of planned .3. Finite k/outer epsilon remain positive; sampled inner covariance reaches -.09983101642. No completed pulse passage or long-time claim. |
| [Initialization build 001](initialization-build-001.json), [002](initialization-build-002.json) and [runtime](initialization-build-runtime.json) | First fingerprint rejected the missing auxiliary transition probe; the complete second target build passes with 100 runtime images and 3,920 numerical inputs at ed34052. |
| [Case-generator CLI](casegen-validation.json) | The documented Regularized command, explicit CubicScaled/C2 controls and default-source Original control emit the requested selectors. No numerical qualification is inferred from generation. |
| [Capped native](capped-native-validation.json), [protocol](capped-native-preregistration.json), [runner](capped-native-runner.json) and [index](capped-native-curation-index.json) | All six severe-depletion cases pass on corrected ed34052: independent complete-ray source/trace/energy audits, two free-k difference levels, and one actual native stage with residual 2.45e-11. This is source-mechanism verification, not general high-order positivity. |
| [Earlier regression JUnit](legacy-ctest.xml), [log](legacy-regression.log) and [campaign runner](campaign-runner.json) | All 33 earlier tests pass in 3,689.71 seconds, CTest parallel12. The enclosing original campaign remains failed because of the separately preserved Path receipt errors. |
| [Corrected-binary Original controls](original-controls-validation.json), [JUnit](original-controls-ctest.xml), [before](original-controls-runtime-before.json) and [after](original-controls-runtime-after.json) | Default Original flat Functions MMS and fresh Neutral evolution pass. Exact-field/traction verification and supported accounting are distinguished from the unsupported full manufactured-forcing energy audit. |
| [Resource accounting 001](resource-adjustment-001.json) and [002](resource-adjustment-002.json) | Explicit numerical/thread/build allocations; conservative combined maximum 38 workers, hard limit 72. |

Both native receipts explicitly retain failed finite covariance gates and
set `candidate_qualified=false`. The dt study holds the physical model fixed;
p and mesh are varied separately. MPI and restart compare every field in
physical mass norms, and restart also compares every native history slot.
The largest native field discrepancy is below 4.1e-14.

The original AA native/MMS attempt stopped on a Python Path serialization
error at checkpoint registration. Successful exit-zero solves are reused
only after regenerated XML and report bytes, solver and actual loader images,
rank counts and the prior producing runtime pin agree. Every independent
audit is rerun under source 628; failed cases and checkpoints are never
replayed. The completion performs the missing checkpoint/restart and
additional MMS solves afresh. Original failed receipts/logs stay on Linux;
the final campaign receipt records their exit codes instead of relabeling
the original campaign as passed.

The [static source-epoch audit](source-epochs.json), generated by
[source-epochs.py](source-epochs.py), records seven immutable Git commits,
11 original frozen files and four unchanged kernel/probe inputs. All four
kernel/probe files are identical from 07363199 through ed34052. Among C++
and header inputs, AA→628 changes only the separate CurvedWallProbe;
628→ed changes only the native startup/signature implementation. This is a
source identity audit, not a substitute for producing binary fingerprints.

The curated runtime summaries retain every image hash, resolved path, actual
loader binding, environment, package/tool identity, mandatory original and
comparison result. Only full `source_inputs`, documentation-tool maps and raw
`ldd` text are omitted; each summary supplies the original file hash/path,
input count and canonical source-manifest hash. Full originals remain on
Linux. The separate mapped probe retains its actual ldd text above.

The capped fixture constrains the actual projected trace coefficients and
their dependence on free local k. Its coarse face-ramp projection has a
pointwise analytic-trace discrepancy .0184441, which remains recorded. No
pointwise analytic equality, unsampled polynomial positivity or broad model
qualification is inferred from the exact coefficient constraint.

Historical [pilot 001](pilot-001-runner.json), [pilot 002](pilot-002-runner.json),
the successful [pilot 003](pilot-003-runner.json), the failed
[Original import attempt](original-controls-attempt-001.json) and
[unsupported energy-oracle request](original-controls-attempt-002.json)
remain separate records. The Original [completion recipe](complete-original-controls.py)
reuses only the successful hashed MMS solve and runs the missing Neutral
control afresh. Historical runner files are retained as evidence of executed
attempts; their existence does not label a failed attempt as passed.
