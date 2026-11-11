# Native ownership correction and requalification

Starting revision is `5a1bc150b72df96bafb2f3cc8b75479064fd4bff`, with a clean
implementation worktree and matching remote. The unrelated main Mac checkout,
Linux primary checkout, CLS/WJ campaigns and every historical evidence directory
remain unchanged. Specifications retain hashes
`056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64`
and `17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40`.

The new authoritative Linux root is
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111`.
Source, build, dependency, frozen runtime, tmp and trajectories are separate.
Numerical libraries use one thread, native execution is serial flat 2D, and
the owner-aware resource guard reserves startup headroom within 72 workers.

## Contract before changes and evolution

The legacy operator is labelled `LegacyCoordinateToleranceV1` in this report;
its checkpoints did not carry a separate ownership tag. The corrected numerical
operator will be `GeneratedVolumeOwnerV1`, separate from the physical resolved
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

Pending corrected production implementation and staged Linux qualification.
