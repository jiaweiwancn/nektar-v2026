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
