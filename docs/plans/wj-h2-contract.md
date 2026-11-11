# WJ-H2: open rough-flat boundary and precursor contract

Reviewed/fetched baseline117e2e6b6dce6f9887b59a17b68c13d081c5ae30; Mac WJ
worktree remains on codex/wj-earsm-awf. Preserve all H1 and WJ7 sources/evidence:
WJ7 RESOURCE_LIMIT/T3 NOT_RUN is closed. Linux root is independently allocated
`/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111`. Native changes use a fresh
source/build, source-built SpatialDomains and verified loader environment.
No installed runtime or other worktree/campaign is changed.

Specifications retain hashes1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124
and0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c.
Keep the literal quad/\\quad erratum and coefficient(d+a*h)/a^2. Original
Cao–Tamura PDF hashadffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df.

## Physical inputs and frozen completion design

Only rough/rough and rough-flat: H_h=1,L_h=2.5,z0=.005, assumed d0=.0875,
delta~6.25,Ztop=20,u_tau=.051,nu=.051*.005/3.97. Distinguish cubes, sand grain,
displacement and fitted aerodynamic intercept. Extract only hollow rough-case
circles in Fig3(a-d), p7, retaining PDF marker centers and labeled-axis
calibration, panel-specific physical height origin and reading resolution.
Visually check overlays; no fabricated probe-error bars or figure profiles.

The virtual crest x=0 is identified approximately with the imposed inlet plane
for this first precursor; domain[0,8] with stations0,2,4,6,8 in H_h units.
The zero-fetch data plane is prescribed, not an independent developed-profile
match. Missing upstream fetch remains a boundary approximation, quantified by
downstream changes. No hill speed-up/reattachment fit is permitted.

Complete measured spanwise variance with weights0,.5,1 between measured vertical
and longitudinal variances. K is half the sum of all three variances. Epsilon
is unmeasured: use a stated finite-depth mixing-length family, calibrated at H
to A(K)/H, with positive free-stream TI=.003 and length-scale sensitivity.
Neither measured k nor epsilon is called exact WJ equilibrium. Evaluate full
WJ covariance and rank-one principal forms across completed profiles/blends
before expensive evolution; retain every inadmissible completion witness.
An alternative must be justified by missing physical inputs, not a q cutoff.

Select the smallest candidate H_w in{.375,.5,.75} containing at least three
digitized overlap mean points above the assumed effective wall. Fit one h_s
using only those rough-flat means and reported wall stress u_tau^2, with K(H)
from the chosen completion. Preserve residuals and check other overlap heights,
h_s<H, finite-nu cutoff, reference-column positivity and molecular fraction.
Freeze the rule before any future hill comparison. H1's values are starting
candidates only. No simultaneous exact matching of measured stresses and WJ
neutral rms1.759/1.192 versus2.15/1.13 is required.

## Boundary equations and implementation scope

New explicit serial flat composite-tunnel WJ Budget mode OpenTunnelV1; existing
prescribed-boundary verification remains the default. Use the normal
IncNavierStokesSolver session/Functions entrypoint, physical coordinates. Q3
differentials, regional continuous Q1 pressures, shared U/k, outer-only physical
epsilon and current dependent H trace remain. No inner epsilon PDE/extra H flux.

For already-integrated conservative volume advection, add exterior fluxes
U*(U.n)-T.n, k*(U.n)-Dk*grad(k).n, epsilon*(U.n)-De*grad(epsilon).n.
Use complete regional physical stress; pi=p/rho+2k/3. On the outlet prescribe
the full mean traction vector, not both pressure and velocity. This removes
the pressure constant nullspace: remove the gauge unknown/row (retaining an
empty report block), not a second pressure-level condition. Physical precursor
uses zero mean traction, scalar zero diffusive flux and explicitly stops on
material outlet backflow; no positive-velocity clamp. This outlet approximation
can cause development and needs a buffer/sensitivity in a future hill case.

Ground: no-slip and modeled zero k flux. Top: impermeable normal velocity,
specified tangential traction (zero for the physical slip approximation),
specified outward scalar diffusive fluxes (zero physically). No top-Dirichlet
mass-flux repair. Inlet differential values are prescribed, ground no-slip
takes precedence, H epsilon belongs to the dependent trace; reject conflicting
independent epsilon corner values. Keep both H-sided outlet owners.

Residual, frozen Oseen/scalar blocks, affine pressure action and independent
budgets must use the same exterior terms. Natural outlet traction contains no
pressure unknown; constant pressure therefore has a nonzero free-velocity
action. Preserve S=D*(alpha*M)^(-1)*G-C with consistent M and roughness G.
One U, pressure, k and epsilon solve and one logical bound per completed attempt;
hard predicates/StopOnBound retain candidates before history commit.

Analytic control is separate: K=qe*u_ref^2, A=u_ref^3/kappa, integrated Budget
inner k, U_H=u_ref^2*J0, physical r=z-d0 in outer log/epsilon (not outer rho),
pi_o=0, pi_w=-K/Ne. Supply the complete reference outlet traction, including
transverse+u_ref^2, and top tangential shear and nonzero epsilon flux. No
experimental z0 is separately imposed on that control. Tolerances account for
nonpolynomial projection. Do not reuse the incompatible H1 vortex as a clean
traction test; any retained curved projection study is passive only.

## Model, reuse and focused checks

Full3D WJ/cubic-N, Original/DevStrainRawWorkV1, raw production and independent
Cd transport remain at c1=1.8,c2=5/9,Ce1=1.44,Ce2=1.82,sigma_k=1,
sigma_epsilon=1.3,Cd=572/6561,qe=81/(2sqrt143),kappa=.381917863546733.
Positive wall nu, Linear full-strain inner stress and uncapped P_A-E_c persist.
No forcing, clipping, regularization, equilibrium restoration or source cap.
Default3000-state guard remains. EP, general3D/MPI, sparse rewrite and WJ7
continuation are excluded.

Qualify any ExactStepV1 extension explicitly on this route: immutable map/mesh
rows only, exact full-state keys, proven scalar-block reuse, time/data generation
and boundary identity invalidation. Keep outlet backflow checks outside scalar-
only cache reuse. Compare Off/On on identical admissible and cutoff/covariance/
principal-failing states, changed boundary data and time. Do not drop inventories
or change tolerances. Profile actual domain.compute/assembly/solve/I/O costs.

Tests: independent physical boundary measures/signs/loads, all H-sided outlet
segments, frozen same-state blocks, pressure columns versus directional
differences, absent redundant gauge, mass/scalar budgets and reason-specific
backflow/corner failures. Re-run relevant small Budget/legacy controls under the
verified source-built library; add boundary/cache restart/rollback checks.

## Finite execution and comparison design (before launches)

Linux one rank/numerical thread, initially one compiler, explicit six library
thread controls, -ffp-contract=off -fno-fast-math, RANS_BUILD_SOLVER=ON. Current
host~251GiB/239GiB available,1.43TB free and22 observed scientific workers; inspect
again before launch and keep all concurrent workers<=72. Build cap7200s;
numerical-study cap21600s including observers/I/O/waits/failures and600s closure;
<=2200 scheduled governing/control units, RSS8GiB, invocation output8GiB,
root32GiB, native invocation<=1500s. Preserve unknown reservations and failed
attempts in a new H2 journal; never touch old ledgers.

Pilot: four fresh physical steps at dt=.02, T=.08, followed by compatible fresh
half-dt eight steps to the same time. These count toward any admitted full run.
Use one horizontal element initially, fitted H/cutoff/roughness-aware normal
knots (target24 cells), confirm actual dimensions. One targeted normal refinement
(target36 cells), same inlet/H/h_s, remains below3000. If a demonstrated
implementation error requires a changed design, retain it before a new launch.

Freeze measured remaining cost with25% headroom and closure before continuing.
Target physical interval T=8, one nominal U_infinity flow-through of L_fetch=8:
400/800 total steps at .02/.01, plus400 refined and up to400 one missing-input
companion after the base. Maximum planned2000+bounded controls. Native blocks
<=32 steps, observations at least every T=.5 and each final common time.
Any hard loss/correction stops immediately. If the cost cannot fit, retain
RESOURCE_LIMIT; do not silently shorten the frozen target. A temporal pilot
failure can motivate one declared dt=.01/.005 design before the full run,
subject to the same step/time caps; it is not an unlimited halving search.

Before inspecting differences, numerical comparison limits are normalized
physical L2 U/k/epsilon2%, wall-stress5%, and matching/cutoff2%; sampled profile
maxima5% (absolute1e-6 near zero), pressure aligned only if a gauge exists.
Use identical physical points/weights and amplitudes U_inf, K(H), epsilon_H,
u_tau^2. A paired/two-level result is sensitivity, not an asymptotic order.
Experimental discrepancy and missing-input spread remain separate.

Development criteria: normalized uncancelled momentum/continuity/k/full outer-
epsilon residuals<=1e-5 and station profile changes<=1e-4 over the last three
equal windows. Otherwise label the attained interval transient. Report T and
T/L_fetch, never infer stationarity from a small stage defect. Preserve all
unreached targets and terminal PASS/OBSTRUCTED/UNRESOLVED/RESOURCE_LIMIT/NOT_RUN
components. No hill prediction or experimental validation claim belongs to H2.

The first production CP0 attempts (no governing steps) exposed and retain input
parser/map-selector errors and an affine native GetCoords endpoint just outside
the certified computational ceiling. The fix canonicalizes only out-of-element
projection-coordinate roundoff within16 binary64 eps times the element scale;
actual out-of-element coordinates and all physical map/cutoff/state guards remain
strict. No interior sample or field is clipped. Source-built Budget and legacy
controls qualify the shared projection change. Measurement packets can request
physical coordinates directly; production code contains no experimental point
list. Preserve producer/observer identities across this development correction.

The first flowing study is closed OBSTRUCTED: coarse four steps at dt=.02 and
fine eight at dt=.01 reached T=.08, with exact Off/On trajectories. The predefined
1966-state normal refinement accepted two steps to T=.04; its T=.06 candidate
failed the in-plane principal predicate with positive covariance and no numerical
correction. No T=8 continuation is admitted. Preserve its complete failed state,
incoming history and accepted predecessor. Authorize one diagnostic refined
half-timestep companion: fresh compatible CP0, dt=.01, at most eight steps to
T=.08, stopping at the first hard/correction exit. Reserve300s and8 units plus
zero-step startup/passive time in the SAME H2 ledger. This diagnoses timestep
sensitivity, does not replace the failed study or admit a new campaign.

New versioned passive observations may select earlier accepted native checkpoints
from their full original producing report, or an accepted predecessor retained
in a failed report under an explicit diagnostic role. Bind exact selected state,
step/time, fixed reference and actual per-entry callback/audit. Terminal native
histories retain full report equality checks; earlier native checkpoint histories
retain their checksum/signature and primary-state checks and are reported as
such. Never use a zero-step report to fabricate a preceding callback or promote
the rejected candidate to an accepted checkpoint.
