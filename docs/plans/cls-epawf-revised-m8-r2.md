# M8-R2: restricted native inner law and admissible implicit stage

Status: in progress. M8 remains **open**. This report will record the actual
bounded stage result, including unsuccessful verification attempts; it does
not authorize terrain calibration or a long evolution campaign.

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
declared numerical thread settings equal to one. With at most three analysis
workers, the allocation is at most 19 workers, below the hard limit of 72.
Producing source, programs, loaded libraries and inputs are fingerprinted.

The two permitted initial changes refine only the outer mesh (6 then, only
if needed, 9 outer cells). The inner mesh and physical inputs stay fixed.
Their preflight counts are 2165 and 2759 unknowns, respectively; the dense
3000-unknown guard remains enforced. Counts are planning data until checked
against native field maps. Both base and higher quadrature are checked.

Only after operator/action verification and initial preflight, the matched
old/new stress problems use the identical initialization, time `.0025`,
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
initialization. Linear has reached its external wall limit; SaturatedR2 is
still running. The Linear target alpha400 was observed in the actual log,
with rounded residuals `62.5985,15.4194,13.1962,13.185,13.1849`, before
interruption during its next Jacobian. No completed target record or final
case report was published. Its target convergence, final scalar state and
physical/accounting gates remain unavailable, rather than false.

The completed Linear alpha3200 auxiliary search stagnated at
`5.334240445678192e-12`, above the unchanged `2e-12` coupled tolerance.
A focused independent terminal replay gives `6.091598802446764e-12` and
also fails that tolerance. Its finite scalar inventory, covariance and fixed
outer-box checks pass: the minimum dimensionless scalar margin is `.675284`,
covariance minima are `.0746028/.176920`, and the maximum normalized outer
gradient departure is `.000532163`. Lambda and b are one, with no active
limiter/cap source volume. These physical-range observations do not turn an
auxiliary residual minimizer into a converged stage or a physical time sample.

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

## Decision gates still open

The terminal R2 decision must distinguish constitutive feasibility, native
weak-operator correctness, convergence to the unchanged coupled residual,
scalar/covariance/outer-box checks on that root, and fully discrete or
physical qualification. A finite-set pass is not positivity everywhere.
The old broad fixed-gamma and invariant-response counterexamples remain
failed; improved numerical convergence cannot reverse them.

Genuine curved MMS for an adopted law, active-cap curved/variable-H coupling,
outer and fully discrete turbulence admissibility, completed perturbation
and long-time refinement, scalable feasibility beyond the dense backend,
and paired developed floor/hill validation and reserved physical scores
remain open.
