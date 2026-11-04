# M8-R2: restricted native inner law and admissible implicit stage

Status: final evidence collection. **No admissible coupled target-stage root
has been obtained. M8 remains open.** The opt-in inner law and scalar-domain
policy are implemented, with scoped native operator verification and an
initial representation inside the fixed finite-set contract. The matched
target searches exhausted their execution budgets without publishing complete
target states. A separately permitted frozen-scalar target search ended in
residual stagnation near the scalar boundary. An independently observed
non-pressure Jacobian-action defect also remains above its unchanged gate.
These results establish neither root absence nor physical qualification.

## Contract and provenance

The reviewed baseline is `5186d63bd825aadd198488d93cff02fc9496a232`.
The Mac branch started exactly at that revision with no tracked changes.
Untracked user literature and HIT material is preserved. R1 and earlier
evidence is immutable. Source and Git work use the Mac checkout; all numerical
results use the designated Linux server over Tailscale. The server's main
checkout remains on its unrelated `codex/hit-inflow-generator` branch at
`f4b70baec10016dd098005bc5a5507b87e4a57bf`.

The [corrected contract v2](cls-epawf-revised-evidence/m8-r2/candidate-contract-v2.md)
fixes the distinct molecular viscosities and dimensional gamma. It retains
the R1 outer certificate and all coefficients. Both implemented pressures
remain reduced pressures; no mechanical-pressure continuity or outer
molecular diffusion is introduced. The [operational policy](cls-epawf-revised-evidence/m8-r2/operational-policy-v1.md)
declares the finite evaluation inventory and the stricter serial, flat,
full-dense prototype restrictions.

The new `RANSWallInnerStress=SaturatedR2` and
`RANSWallNonlinearPolicy=ScalarFeasibleR2` selectors are opt-in. Linear stress
and the Legacy solve remain defaults. The selected stress is shared by
momentum, production (including weak divergence), traction, energy and
covariance. Scalar diffusivities, analytical auxiliary profiles, capped
budgets and the dependent epsilon trace retain their specified definitions.
See the [caller audit](cls-epawf-revised-evidence/m8-r2/inner-stress-caller-audit.md).

## Bounded execution

Heavy source/build epochs, fields, logs, matrices and raw reports live below
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r2`, outside the checkout.
Compilation uses 16 workers. Numerical programs use one MPI rank and all six
declared numerical thread settings equal to one. Native searches are sequential;
concurrent independent analyses receive explicit individual CPU allocations.
The recorded allocation and process checks remain below the hard limit of 72,
including the idle MPI management threads observed in native processes.
Producing source, programs, loaded libraries and inputs are fingerprinted.

One targeted outer refinement was sufficient: outer6 has **2165 actual
unknowns**, verified against native field maps. Its dense matrix occupies
37,497,800 bytes. The second permitted representation, outer9, was not run;
its prospective count is 2759. The inner mesh and physical inputs stay fixed,
and the dense 3000-unknown guard remains enforced. Base and higher quadrature
are independently checked.

After the scoped inner-operator, pressure-action and initial-state checks, the
matched old/new stress problems use the identical initialization, time `.0025`,
history, forcing and tolerance `2e-12`. Auxiliary mass coefficients are
`3200,1600,800`, followed by the exact target `400`. Auxiliary solves never
advance history. After an auxiliary failure, remaining auxiliaries are
skipped and the target is still attempted once. The maximum is 12 Newton
updates per level. No additional alpha search is authorized.
Each native case has an external 14,400-second wall limit. A wall-limit exit
does not supply a nonlinear terminal classification or a missing target vector.

Two short native runs are conditional on a new-law target passing every
required gate. Otherwise, if both targets remain unresolved, at most one
frozen-velocity/pressure scalar-stage diagnostic is permitted. A failed
search does not prove absence of an admissible root.
That diagnostic freezes represented velocity and the native **reduced
pressure** coefficients Pi. Mechanical pressure continues to satisfy
`p/rho_f=Pi-2k/3`, so all k-dependent mechanical-pressure, sensor, source,
normalization and projected-trace terms remain in its scalar residual.
It solves a declared subsystem and supplies no coupled-stage or history
advance.

## Retained unsuccessful attempts

The first isolated build failed a strict indentation warning in a kernel
verification probe. The second failed the same warning class in a native
fixture. The third build at `1ee3a8454c857c807bcb6911024186afb9d22721`
compiled the full solver and probes successfully; compiler flags and equations
were not relaxed. All build logs remain available.

The first inner verification attempt completed the compiled-kernel checks,
then stopped at a full-interface-traction discrepancy of
`1.891880318360961e-5` against `3e-10`. Its original raw output and failed
receipt are retained. Diagnosis and any correction are recorded separately.

The first pressure verification produced both native reports, then failed
to serialize a NumPy boolean in its Python report. This is not a scientific
pass. A separate analysis-only replay checks the exact retained reports and
all independent derivative gates, without silently rerunning native cases.

The legal depleted operator fixture subsequently exposed an overflow in the
historical Python oracle's nested exponential. R2 now uses an isolated stable
indexed evaluation, preserving every historical module and the original
production C(q). No q or field is clipped. The positive-real-axis omitted
coefficient/tangent bounds are below `3e-64` and `4e-62`, respectively. Saved
native outputs were reused. The exact-zero fixture also exposed conditioning
of its diagnostic maximum a/c at local k about `1.2e-13`; its original
coefficient-derived comparison remains **unresolved at 3e-10**. Independent
actual-primitive/kernel, stress, covariance, weak-load and energy checks keep
their original thresholds. The [inner result](cls-epawf-revised-evidence/m8-r2/inner-verification-result.md)
records this limitation and all retained failures.

## Completed operator and initial-state evidence

Producing source `9f20b3053c303e8805f5183f668a8ec468d9656b` builds the
complete Release solver and probes. The current branch's subsequent changes
are independent helpers, their diagnoses and evidence; the producing C++
source remains unchanged. Three focused Legacy regressions pass. The scoped
inner replay completes 4366 checks over 38 compiled cases and eight native
fixtures, with zero native reexecution in its final replay. At 21 exact-zero
points the selected R2 covariance is exactly zero, while the matched Linear
control remains negative. A retained cap-active fixture independently verifies
the full weak residual, exact mass, constraints and projected matching trace.
See the [hash-bound inner summary](cls-epawf-revised-evidence/m8-r2/inner-verification-result.json).

The [pressure localization](cls-epawf-revised-evidence/m8-r2/pressure-action-origin.md)
resolves the independent action discrepancy through exact affine assembly
and complete analytic sensor/source/trace contributions. The absolute
physical-norm action errors are at most `9.79e-11`, below the unchanged
`1e-7` gate. Historical cached/fresh errors and finite-difference cancellation
and truncation failures remain recorded. This is derivative evidence, not a
claim that pressure differencing caused or repaired the old transition failure.

The first outer6 native initialization succeeded, but its independent helper
selected the inner geometry for one outer epsilon evaluation at H. The
helper side selection was corrected; the original report was reused in place
and the failed audit preserved. Three further initializations completed the
four-case matrix. All four complete initial vectors are bitwise identical,
SHA256 `f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca`.
The 3540 independent checks and full runtime/source/input provenance pass.

| Initial finite-set property | Linear | SaturatedR2 |
|---|---:|---:|
| Actual unknowns / source quadrature pairs | 2165 / 8,16 and 12,24 | same |
| Maximum outer gradient departure (limit .001) | .000500194924459 | same |
| Minimum outer k and epsilon box margins | .0100 / .0032524731 | same |
| Minimum inner / outer covariance | .0772933168271 / .176807110004 | .110802957528 / .176807110004 |
| Full inventory sampled K/K0 | about 1 | same |
| Normalization lambda / limiter / cap | 1 / inactive / inactive | same |
| Maximum whole-inventory inner a/c | 1.09303567 | same represented inputs |
| Saturation-active source volume, base / high (inner volume12) | not selected | .04210427 / .04312037 |

Outer6 covers the unchanged fixed outer certificate at both quadrature
levels on the declared finite sets, including the actual and projected
dependent trace. **Outer9 is not needed or run.** The stricter operational
sampled-K guard has margin about .01. This is initial coverage, not
stationarity or a positivity theorem. The projected field activates
saturation although the continuous analytical reference is inactive.

The starting R1 improved mesh has 1571 unknowns (19,744,328 dense-matrix
bytes), with the same inner intervals and only three outer intervals.
Its retained one-sided outer witness at r=8/3 has normalized gradient-entry
departure `.003412338694780`, exceeding `.001`; covariance and positive
scalars alone therefore do not cover the fixed certificate. Its raw
streamwise/k norms are about `9.86323/62.5985`. These are reused R1
representation observations, not a new R2 trajectory or a replay of the
different p7 zero-step mesh. Outer6 reduces the measured gradient departure
to `.000500194924459` while increasing the actual unknown count by 594.
At base quadrature the selected inner stress changes integrated mean stress
work and production equally by `-.0123957428841`, with zero exchange defect.
The unanchored velocity/k norms remain nonzero: about `9.86324/62.59848`
for Linear and `8.22728/62.60960` for R2. The [preflight notes](cls-epawf-revised-evidence/m8-r2/preflight-notes.md)
and [machine table](cls-epawf-revised-evidence/m8-r2/preflight-002-summary.json)
retain gradients, scalar/trace margins, exact quadratures, source accounting,
quadrature sensitivity and costs. Native initial ScalarDomain calls take
about .19s at base quadrature, with the raw residual/constraints/norms about
.13s; these do not measure an entire Jacobian.

An independent adapter precheck matches all 40,627 native fixed-inventory
scalar observations (39,826 unique descriptors) and independently assembles
the polynomial weak-continuity/gauge block. Its direct residual error is
`2.11e-15` in the physical norm, below the unchanged `1e-7` threshold. The
[retained receipt](cls-epawf-revised-evidence/m8-r2/analyzer-static-result.json)
contains 500 checks and no native calls. This verifies the independent field
map and affine block; current stage Jacobians still require their own checks.

The matched actual startup-stage searches use this common admissible
initialization. Both cases reached their external wall limits, with producing
source/runtime fingerprints unchanged. Both actual alpha400 targets were
observed in their native logs. The Linear target has rounded residuals
`62.5985,15.4194,13.1962,13.185,13.1849`, before
interruption during its next Jacobian. No completed target record or final
case report was published for either case. Target convergence, final scalar
state and physical/accounting gates remain unavailable, rather than false.
SaturatedR2 began its exact target after its failed auxiliary search; the
case limit interrupted its first Jacobian, whose last logged column was 256. Its initial
target residual was `62.6096`. A fresh-Jacobian build is not a completed
nonlinear update or a proof that the unchanged target has no admissible root.

The [matched native summary](cls-epawf-revised-evidence/m8-r2/matched-stage-summary-001.json)
retains every complete solve, trial fraction, block norm, scalar margin,
matrix hash and wall-limit command. Missing final-report flags stay null;
the separate strict assessment establishes actual target attempts from the
logs and checks that no complete target-convergence evidence was discarded.

The completed Linear alpha3200 auxiliary search stagnated at
`5.334240445678192e-12`, above the unchanged `2e-12` coupled tolerance.
A focused independent terminal replay gives `6.091598802446764e-12` and
also fails that tolerance. Its finite scalar inventory, covariance and fixed
outer-box checks pass: the minimum dimensionless scalar margin is `.675284`,
covariance minima are `.0746028/.176920`, and the maximum normalized outer
gradient departure is `.000532163`. Lambda and b are one, with no active
limiter/cap source volume. These physical-range observations do not turn an
auxiliary residual minimizer into a converged stage or a physical time sample.
The [compact terminal record](cls-epawf-revised-evidence/m8-r2/terminal-linear-auxiliary-result.json)
binds the independent reconstruction and actual producing source/runtime,
including full traction, k flux, separate weak/strong continuity, constraints
and energy/source accounting. These are independent terminal observations;
the absent final case report supplies no native postqualification table.

The SaturatedR2 alpha3200 search also ended in
`feasible_residual_stagnation`, after 12 accepted internal updates. Its
independent [terminal replay](cls-epawf-revised-evidence/m8-r2/terminal-saturated-auxiliary-result.json)
completed in 229 seconds. It passes the same finite physical-range checks,
while its coupled residual remains above tolerance. These two auxiliary
observations separate scalar/covariance/applicability from convergence:

| Completed auxiliary, alpha3200 | Linear | SaturatedR2 |
| --- | ---: | ---: |
| Actual / independent coupled maximum | `5.33424e-12 / 6.09160e-12` | `4.59464e-12 / 5.47327e-12` |
| Actual coupled tolerance met, `2e-12` | no | no |
| Actual / independent scalar-block maximum | `6.07507e-13 / 7.98783e-13` | `6.70615e-13 / 7.98319e-13` |
| Minimum dimensionless scalar margin | `.675284` | `.675182` |
| Inner / outer covariance minimum | `.0746028 / .176920` | `.0820016 / .176920` |
| Maximum normalized outer gradient departure | `.000532163` | `.000532248` |
| Actual/target outer epsilon/epsilon0 range | `[.9832624,1.0000206]` | `[.9832624,1.0000206]` |
| Common finite box, covariance and scalar conditions | pass | pass |
| Minimum b / lambda; limiter/cap source measure | `1 / 1; 0 / 0` | `1 / 1; 0 / 0` |
| Saturation-active source measure, of 12 | `0` | `.0456635` |
| Complete internal updates / valid non-descent trials | `8 / 33` | `12 / 51` |

The recorded auxiliary trials contain no scalar-domain rejections or residual
evaluation exceptions; none of the fraction-to-boundary limits is active.
Late dimensionless steps are about `6.6e-14` / `4.7e-14`; linear residuals
are about `1e-26`. These observations classify valid non-descent and
stagnation in a fresh-Jacobian policy. They do not establish a roundoff
floor, a causal explanation involving the inaccurate k action, or a linear
factorization failure. Neither failed auxiliary candidate supplies the next
target guess; alpha400 starts from the source-bound h0 when no auxiliary
level converged. Native history remains unadvanced by this probe protocol.

The full auxiliary-path audit finished but failed five moving-ray lookup
checks. A separately bound terminal diagnosis found a one-ULP radius
difference crossing a 12-decimal key boundary. The [association correction](cls-epawf-revised-evidence/m8-r2/ray-coordinate-association.md)
preserves the original coordinate/weight tolerances and adds a stricter raw
coordinate safeguard when keys differ. Eight focused Linux fixtures and the
217-second terminal-only replay pass. The original failed full-path audit
remains failed and has not been rerun. The producing C++ and binary are
unchanged; this analyzer indexing issue is not a cause of native stagnation.

The retained first auxiliary Jacobian exposes a separate numerical action
failure: its outer-k direction differs from an independently assembled
analytic reference by `7.15997e-7`, above the unchanged absolute `1e-7`
gate. The largest contribution lies in outer momentum load rows, with a
smaller roughness-force contribution. Pressure/gauge and affine continuity
checks pass. Two complex steps and both invariant branches agree within
`8.40e-15`; finer whole-residual differences also show increasing subtraction
error. The [initial-action diagnosis](cls-epawf-revised-evidence/m8-r2/initial-j0-action-diagnosis.md)
preserves the original failed gate, independent controls and exact source/
matrix hashes. This initial numerical derivative defect is separate from
constitutive feasibility and does not establish the cause of nonlinear
stagnation. The bounded native comparison remains diagnostic; physical
qualification requires all independent action gates.

The same analytic defect is bound to the Linear target's actual initial
alpha400 matrix through an independently checked exact mass change. The
SaturatedR2 initial auxiliary matrix separately fails the absolute
finite-difference comparison, with errors `7.35528e-7`, `7.86320e-7` and
`1.64307e-6` on the three-step ladder. That is not a second analytic
component diagnosis. Both formal outer invariant branches and the active
sensor/trace fixtures remain part of their respective checks. Neither the
small relative errors nor accurate linear solves replace the absolute gate.

## Frozen-scalar localization and matched decision

The one permitted frozen-scalar diagnostic used the same SaturatedR2 h0,
actual target alpha400, time `.0025`, unchanged tolerance and zero forcing.
It held represented velocity and native reduced pressure Pi fixed, retaining
the actual scalar residual and all scalar-dependent sensor, normalization,
lift and projected-trace constraints. It never advanced native time/history.
The complete record terminates as `feasible_residual_stagnation`.
Its six recorded accepted updates have fractions
`.681356272, .100605631, .000594326, 2.97374e-6, 1.48688e-8, 7.43449e-11`.
Each is limited by the fixed-inventory fraction-to-boundary estimate; all
six trials are evaluable. The limiting native-GLL location is
`(1.6046496089512394,0,0)`. Fresh Jacobians are used at each update.
The native level's `iterations=5` is its final zero-based iteration index;
the complete JSONL record contains six updates.

The native scalar maximum remains `15.545717638361056`; independent terminal
reconstruction gives `15.545717638361037`. The k block dominates, while the
epsilon block is `.06176892447`, the unsolved streamwise block is
`8.127147237`, and the dependent-constraint norm is `1.09764e-9`.
The scalar subsystem alone therefore stagnates near its evaluation-domain
boundary under this policy. This localizes a difficulty beyond just searching
the velocity/pressure unknowns; it does not distinguish a reaction,
diffusion or nonlocal-trace cause, prove root absence, or remove the known
action-accuracy limitation.

The [independent terminal packet](cls-epawf-revised-evidence/m8-r2/terminal-frozen-scalar-result.json)
binds the actual complete record, checkpoint, history, native audit table,
producer and analysis source. Its eleven phases completed in 226.5 seconds,
with one independent residual reconstruction and no native rerun or current
Jacobian/path analysis. The minimum inner k is `7.97763e-15` natively and
`8.04319e-15` independently. The latter is below the verifier's
`9.14768e-15` arithmetic-roundoff indicator (`sign_resolved=false`). That
indicator is not a rigorous error bound; the evaluated positive value is
not a certified positive margin or a polynomial-everywhere result.
Evaluated inner/outer covariance minima are `1.57053e-15/.176810977`.
The unchanged finite outer box is covered, with normalized gradient
departure `.000500194924462`, k/K0 in `[.999995974,1.000012069]`, and
epsilon/epsilon0 in `[.983296232,1.000047449]`. These are observations on a
failed diagnostic candidate, not checks on a converged RANS stage.

Lambda remains one. The limiter-active source measure is `.00103271248`
and minimum b is `.832365063`; cap-active measure is zero. Independent
integrated E and D are both `8.43645916235`, with redistribution integral
about `-1.60e-16`. Active saturation measure is `.05983489225` of inner
volume 12; selected mean-work and production changes both equal
`-.11843748537`. Including chi=0, the actual matching trace remains
the projected `lambda*Epre(H)`. Limiter activity is not normalization-cap
activity, and these observations supply no cap-induced epsilon-attenuation
explanation.

| Layer / actual target status | Linear inner law + ScalarFeasibleR2 | SaturatedR2 + ScalarFeasibleR2 | Frozen SaturatedR2 scalar diagnostic |
| --- | --- | --- | --- |
| Inner constitutive result | Retained linear law and historical failures | Compiled radial-law/derivative/zero-input checks pass in declared scope | Same research law |
| Native selected weak operator | Existing control | Active fixtures, full stress/production/traction and capped constrained residual verified; diagnostic maximum a/c limitation retained | Terminal native/independent residual and budget agreement |
| Initial fixed finite-region coverage | Pass, outer6 / 2165 unknowns | Pass, identical h0 | Same h0 |
| Pressure/gauge and affine continuity actions | Pass unchanged gates | Pass unchanged gates | Frozen rows are not full coupled actions |
| Other current-action qualification | Fail: independent analytic k-action error `7.15997e-7 > 1e-7`, including target initial matrix | Fail: initial auxiliary FD comparison; no complete target-matrix qualification | No new independent current-action/path check |
| Actual alpha400 coupled residual converged | Unknown: wall limit; no complete target record | Unknown: wall limit during first target Jacobian | No: scalar `15.5457`, coupled equations not solved |
| Scalar/covariance/fixed-box checks on a converged coupled target | Unavailable | Unavailable | No converged coupled target; terminal finite evaluations admissible/covered, near-zero sign margin unresolved |
| Complete target source/energy/traction qualification | Unavailable | Unavailable | Independent/native accounting agrees; nonzero residual and interface defects remain |
| Native physical time/history advanced | No advancement by probe protocol; final report absent | Same | Final report confirms time0, steps0, original history retained |
| Fully discrete or physical qualification | Not established | Not established | Diagnostic only |

Accounting agreement must also be distinguished from balance at a root.
For the frozen candidate, independent normal/tangential traction-jump RMS
values are `.001195289/.000377060`, k-flux-jump RMS is `1.40175e-5`, and
integrated reduced-pressure jump is `-.4492470515`. Strong divergence is
zero for the frozen represented reference velocity, while the separate weak
pressure/gauge norms are about `5.05e-16,6.07e-17,1.03e-16`.
The complete energy identity has direct remainder `2.58e-14`, but free weak
work remains `.748357889`; a small accounting identity error does not make
that work zero. No mechanical-pressure continuity condition is imposed.

The frozen supervisor failed with `BrokenPipeError` while printing progress.
The child continued and published a complete native report and solve record.
The failed raw receipt is retained unchanged. A separate post-completion
collection binds those artifacts and a fresh successful source/runtime
fingerprint; lost process exit status and elapsed wall time remain null.
The solver log's `3033.31s` Execute timing is not substituted for a
supervisor-observed process duration. The supervisor's wall-cap enforcement
after its failure is unverified. This collection is not a native retry.

## Recommendation and open gates

Retain SaturatedR2 and ScalarFeasibleR2 as opt-in research prototypes. The
bounded R2 result is **an unresolved admissible coupled-stage-root problem,
with a localized remaining non-pressure action-accuracy defect and expensive
dense searches**. Fix and independently verify that action before another
exact-target search; preserve the residual, tolerance and finite inventory,
and publish internal updates durably so an external stop cannot erase the
target path. Reserve a subsequent bounded solve for the exact alpha400
target rather than consuming its budget on nonconverged auxiliaries.
This is a recommended next task, not an additional R2 campaign.

The five layers remain separate: constitutive feasibility, native operator
correctness, convergence of the true stage residual, physical checks on that
root, and fully discrete/physical validation. A finite-set pass is not
positivity everywhere. No short dt/dt/2 evolutions were launched because the
new-law target did not pass the prerequisite gates.
The old broad fixed-gamma and invariant-response counterexamples remain
failed; improved numerical convergence cannot reverse them.

Genuine curved MMS for an adopted law, active-cap curved/variable-H coupling,
outer and fully discrete turbulence admissibility, completed perturbation
and long-time refinement, scalable feasibility beyond the dense backend,
and paired developed floor/hill validation and reserved physical scores
remain open. The lambda-to-zero / outer high-Re epsilon continuation remains
an unresolved model choice: the prototype supplies no floor or continuation
through epsilon=0. M8 remains **open**.
