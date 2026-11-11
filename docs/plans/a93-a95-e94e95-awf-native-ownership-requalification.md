# Native ownership correction and requalification

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

Fresh staged numerical qualification remains in progress; no corrected-operator
duration or response qualification is claimed by the implementation checks.
