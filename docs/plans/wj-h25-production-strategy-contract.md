# H25: measured production execution strategy

Frozen before new optimization outcomes. Reviewed/fetched branch head is
`a60d6b3ee716ea14d388e8bf6e27dff6eccd63a1`; the WJ worktree is clean and no
H25 root/job was present. H13–H24 remain sealed. H24's mesh/no-build limits were
local to H24; this package permits qualified numerical implementation changes.

The goal is the seven H25 closure criteria, with PRODUCTION_STRATEGY_READY only
when the complete conservative physical-time forecast satisfies
`1.25*T_case <= 86400s`, equation work exceeds50% (target70%), at least95% of
wall cost is attributed, numerical gates pass, and a real pilot/restart/storage
qualification and commits A/B/C are pushed. Otherwise iterate within the finite
development allocation and publish a substantiated TARGET_NOT_MET decision.
Neither a profile alone nor a modest speedup is scientific closure.

## Hardware, source, and finite development allocation

Source/Git: branch-owning Mac `.worktrees/wj-earsm-awf`. Authoritative Linux:
`jay-Super-Server`, two AMD EPYC7763 sockets,128 exposed single-thread CPUs,
Ubuntu GNU11.4/private-loader convention. The documented USB endpoint timed out;
use authenticated Tailscale. The Linux mount alias is currently
`/media/jay/新加卷1/nektar-v2026` (HEAD f4b70baec10016dd098005bc5a5507b87e4a57bf,
unrelated developer-guide modification preserved). Build/run only an isolated
source capsule on healthy ext4 under `/home/jay/build/wj-h25-production-20261116`.
The installed runtime and old evidence are read-only inputs, not new sources.

Development ceilings, not one-case forecasts:14400s total elapsed operations,
2400s build,7200s native/component/benchmark/pilot,3600s passive analysis and
transport. Setup120s and closure600s are fully charged within total; initial
access checks are inside setup. Per operation: build<=900s, native<=900s,
passive<=300s. Unknown/failed operations retain their reservation charge.
At most200 WJ governing attempts and twelve WJ advancing process blocks, including
reference, perturbation, screening and restart work; stock loops have separate
counts and <=120s per segment. No automatic uncertain native relaunch. Admit
remaining work with25% prospective headroom, protect closure and final read.

One serial scientific native process and a bounded observer, six numerical
library thread controls1; build<=8 and total host allocation<=72. Science
aggregate RSS<=8GiB, build<=16GiB. New source/build/runtime<=8GiB, new science
output<=8GiB, free-space reserve>=15GiB before admission. Preserve unrelated jobs;
only owned groups can be terminated. No full505-unit trajectory is launched.

## Scientific invariant and reference comparison

Keep sole6mm primary, d0=3.5mm, interface9.5mm, complete stored roughness fit,
full-height inlet/initial reference, top traction/epsilon flux, all outlet and
interface laws, Original full3D WJ/DevStrainRawWorkV1/raw full-stress production,
Cd, high-Re outer and finite-nu budget AWF, nonconstant transported inner k,
current-K cutoff/dependent epsilon, IMEX1/eta0 and one flow→k/bound→epsilon sweep.
No EP, k reset, forcing, floors, physical regularization, height change or hill
evolution. q_eq, actual qh0, fixed reference/defect and restarted state are distinct.
Horizon505.02466521278404, terminal window[447.8483832973273,505.02466521278404]
stay fixed. Continuation from .04 and complete fresh-start costs are separate.

Unchanged reference W0:700 Q3/Q1 elements,24688 states,14 columns,16+34 normal
intervals, q12/normal16, own-dt qualified H24 state at .04. ExactStepV1 On,
ImmutableRowReuse Off, CheckpointAuditV1/BoundedV1/StopOnBound, all hard guards
and residuals active. Benchmark two successive warm advancing segments, not
repeated identical-state evaluations. Separate startup/factors and final reads.

S0 is stock VelocityCorrectionScheme, no RANS/no SVV, IMEX1, Q3 velocity and its
supported continuous Q2 pressure, positive Kinvis equal to stored molecular nu.
Materialize W0's PHYSICAL mesh (x unchanged, physical normal coordinate from the
composite map), preserving all700 cell knots and9.5mm interface location. Verify
coordinates, determinant/physical measures497.8125 total,3.75 inner,494.0625 outer.
Stock uses nonzero inlet/initial velocity, no-slip ground, fixed top tangential
velocity and normal impermeability, supported natural outlet velocity with
pressure level. Its top/outlet/viscosity/global-pressure differences from WJ are
explicit: it is a cost reference, not a RANS physical substitute. Equal-resource
comparison first uses dt=.0025 and equal physical progress on the same mesh.
Stock setup and warm advances are separated; include repeated segments and CPU,
RSS/output/variability. S1 only after a diagnosed S0 stabilization failure, with
explicit kernel/cutoff/coefficient/reference-speed and separate sensitivity.

An additional common frozen-operator component fixture must use identical basis,
quadrature, matrix coefficients, free constraints/RHS/tolerance. Full WJ frozen
operator timings remain separately labeled. Stock-versus-WJ elapsed differences
alone are not identified as RANS equation overhead.

## Predeclared numerical gates and candidate ordering

Unchanged arithmetic: exact IEEE state/history/row and failure parity. Changed
arithmetic: scaled operator/action/residual errors<=2e-11; state/block trajectory
reconstruction floor4e-11. Native/independent principal replay2e-11. Retain the
existing hard domain/root/scalar/ground/outlet/covariance/transverse gates and
separate normalized principal-quality threshold.002068763923454. Every current
regional Q1 moment<=2e-15; evolved strong divergence remains separate. One bound
call/sweep, zero activations/material changes for a qualified candidate; active
adverse fixtures must preserve first failure, location and correction mathematics.
No tolerance is tuned after seeing a candidate's outcome.

For mesh/quadrature/time differences evaluate BOTH fields on identical physical
Gauss8 and fixed/witness/face-owner support, split at actual faces and physical
wall/cutoff/roughness breakpoints; preserve every owner. Regional physical RMS
and finite-inventory MAX are separately limited: U1e-4 (scale1), k1e-3 (K),
R/pi1e-3 (u_tau^2), G1e-2 (max(||G_eq||,.051/19.9125)), outer epsilon1e-2
(positive local epsilon_eq);4e-11 scaled floor. Preserve signed components,
inner dependent destruction separately, wall traction1e-3 relative to u_tau^2,
current-K/trace/budget and interface flux reports. Same-owner principal action
differences, including each member's minimizing direction, must stay within
`.1*(a_finer-.002068763923454)` with resolved allowance; no pressure regauging.
Weak component validation keeps2e-9+2e-8|term| and regional4e-9+2e-8|sum|.
These numerical budgets are provisional, not wind-tunnel uncertainty.

Implement exact cost removal first on700cells, with dynamic invalidation and
changed-k/epsilon/cutoff/time/owner/adverse-state coverage. Do not re-enable the
slower H17 generic cache. Then screen the recorded parent-knot350-cell Q3/Q1
mesh (14×[8+17]), preserving topology features/physical domain; p reduction only
if real mixed-space support is implemented. Quadrature is a separate candidate,
never confused with degree or the native8×8 representation. Initial trial bounds
are q12/16 baseline, q8/12, then q8/8 only after budget/source checks. Each variant
has its own signature, physical state initialization/transfer and native history;
ordinary wrong-signature restart must remain rejected.

Temporal candidates are .005,.01,.02,.04,.08, screened upward against a smaller-dt
common-time reference, without treating no-crash as accuracy. Freeze exact
initial field/perturbation and common times before each screen. The perturbation
must have finite nonzero velocity/gradient and1% k/epsilon-scale changes, vanish
on prescribed boundaries, preserve required traces and pass native admission;
it must drive changing K/epsilon/cutoff. It is a separately identified transient
fixture, never an equilibrium repair. Representation, h, q and dt errors are
reported separately; no winner is selected from favorable near-equilibrium drift.

## Measurement, operational pilot and publication

Disjoint taxonomy: genuine equation evaluation/assembly/matrix/factor/solve/
required state updates; guards/bound detection and verification; materialization/
allocation/copy; instrumentation; full scientific audit; observer; I/O/recovery;
startup/restart/terminal. Attribute>=95%; resolve the native.implicit exclusive
remainder with extra scoped/native samples. Nested scopes never sum as disjoint
work. Measure instrumentation overhead and warmed advancing variability, not
only cache hits. Profiling permission failure is retained; no host changes.

Candidate production monitoring starts from physical recovery.025/science.25/
terminal.125 scenarios. A selected schedule must be compared to denser auditing
on an evolving perturbed trajectory with explicit event escalation; every retained
scientific checkpoint is freshly qualified. Bound resident support and archive
history; preserve current witness coverage. Compact telemetry and recovery
rotation must have synthetic long-sequence tests plus a real pilot/restart/
strict read. Initial configuration/runtime are pinned once. No exact unsaved
extrema, nonlinear averages or integrated budgets are inferred from sparse data.

Before numerical execution commit this contract. Checkpoint A adds completed
fair stock/W0 benchmarks and disjoint attribution, then continue. Material fixes
are committed before executing their source capsule. Checkpoint B publishes
validated improvements and same-mesh cost/failure parity, then continue. C seals
all candidate decisions, real operational pilot, measured complete forecast,
closure matrix and CLOSED finite campaign ledger, and is normally pushed.

The complete forecast includes startup, accepted advancement, audits/observer,
recovery/output, restart provision and terminal/export; add fresh spatial/native
initialization and0→.04 for a fresh case. Apply25% contingency and separate genuine
equation-work share for ordinary steps AND complete forecast. Hill refinement,
curved-route cost and application budgets remain later requirements. The7mm
observations remain inside the9.5mm matching plane; performance does not validate
that physical comparison. If a concrete unresolved obstruction remains after
bounded iteration, publish TARGET_NOT_MET, the cheapest accuracy-qualified case,
full cost, failed approaches and one precise next implementation change.
