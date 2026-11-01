# CV3 flat-floor calibration protocol

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: computational experiment execution
- Origin Date: 2026-11-01
- Verification Status: bounded design and numerical controls executed; CV3 development/freeze gate unmet
- Version Label: cao-tamura-cv3-protocol-v1

This implements CV3 of the [calibrated-validation plan](cao-tamura-calibrated-validation.md).
The frozen CV1 data and CV2 reference/observation decisions remain in force.
Only flat calibration series 3a/3b/3c and reported flat scalar constraints may
enter selection. Reserved 3d is evaluated only after the selected parameters
are frozen. Hill predictions/targets never enter this campaign.

## Physical and numerical definition

Use separate smooth S and rough R empty floors, the canonical standard and
CLS closures, retained molecular diffusion, a 0.8 m roof and physical normal
H/P/N = 11/5.5/22 mm. Smooth AWF ks is exactly zero. Rough ks is an identified
equivalent sand-grain parameter; cube height, d and z0 remain distinct.
Integrated strip continuity and physical Yap distance are retained.

The primary calibration uses original AWF for both closures. A common set of
inlet/surface parameters is selected across closures for each surface. Transfer
that set to EP-AWF and the compatible standard log wall; do not retune the
inlet to conceal wall-policy differences. The policy-required H..N epsilon
collar is declared separately from the common outer/inlet parameters.
Any closure-conditioned alternative must have its own version and label;
none is implied by choosing a common fit.

For common physical boundary data across closures use the standard reference
Q_in=10/3 and kappa_in=sqrt(0.48*1.3/Q_in) only to parameterize the
unmeasured inlet k/epsilon. Each closure retains its own canonical Q/kappa
and its required H..N epsilon collar. These inlet reference numbers do not
change the closure coefficients.

The inlet mean is initialized from the allowed Fig. 3a log-polynomial, with a
finite-depth outer blend. Calibration controls are incoming thickness,
inferred dissipation scale, incoming k amplitude and rough ks. These change
boundary data, not canonical closure or smooth-wall constants. U_inf/nu uses
the frozen flat U(Hh)/U_inf interpolation and reported Re_H. The nominal target
U(Hh)=1 must be checked at the developed incident station, x=0. Absolute
dimensional wind speed is not inferred from hill amplitudes.

A straight-floor mesh may use fewer streamwise columns than the hill mesh;
all degrees/columns/normal layers and case hashes are recorded. Initial
throughput/CFL pilots select a usable timestep without relaxing physical
residual or positivity gates. Mesh/time convergence, facility sensitivities,
long-time/3D stability and production sampling remain CV4 gates. Reopen a
frozen calibration if those checks change its incident observables materially.

## Declared search and selection

Use a bounded coordinate design rather than claim global optimization:

| Control | Centre | Lower / upper perturbations |
|---|---:|---:|
| Incoming thickness (m) | 0.21 | 0.17 / 0.25 |
| Dissipation multiplier | 1 | 0.5 / 2 |
| Incoming k multiplier | 1 | 0.75 / 1.5 |
| Rough ks (m) | 0.004 | 0.002 / 0.008 |

Run the centre and one-coordinate perturbations, retaining unsuccessful or
inadmissible points. Select a common surface candidate using both closures'
results. If a further bracket is required, append its rationale and bounds
before running it; do not erase the first design or failed results.
Parameter sensitivities/flat directions and any bounded-search limitation
remain part of the result. A best point at a bound is not an identified
interior optimum.

For parameter selection use the predeclared common resolved band z>=N=22 mm.
Keep the frozen height weights, report supported/excluded weight and retain
full-range/conditional profiles separately. Do not turn below-H stress/rms
gaps into full-range scores. Normalize mean, shear, wall friction, nominal
thickness and longitudinal-rms residuals by the planning engineering scales
3% U_inf, 10% u_tau^2, 5% u_tau, 5% nominal delta, and 15% experimental
weighted rms respectively, with equal family weights. These are selection
scales, not measurement standard deviations or statistical confidence limits.
Publish each component; aggregate improvement cannot hide one component's
failure. Missing/negative variance and unavailable support are retained.

Measure friction from native wall traction independently of mean-profile
fitting. Keep the two one-sided numerical traces at x=0 and quantify their
difference. Report incident Re_H, U(Hh), mean/shear/rms errors, delta99 and
delta95, and outer turbulence separately. The source's approximate delta has
no stated edge criterion: delta99 is a declared analysis convention and
delta95 is its sensitivity, not a new source measurement.

## Developed flow and uncertainty

Advance actual native URANS fields with full-history restart and physical
time retained. A prescribed or newly prepared profile cannot count as a
developed precursor. For steady calibration candidates, require at least
three domain flow-through times before assessment and two successive
one-flow-time extensions with small incident/full-field drift. Initial drift
limits are 2e-4 U_inf for mean velocity, 2e-3 u_tau^2 for shear, 1e-3 relative
k/epsilon at every native quadrature point, and 0.2% friction/thickness. Extend a drifting case; a timeout or
growth is an incomplete candidate, not a steady one. Record all stages,
positivity/guard counts, normal/mass budgets and source residuals.

The developed-steady check supports provisional flat calibration on its
declared discretization. It does not replace CV4 stationarity, uncertainty
or a perturbed 3D stability demonstration. Use mean instantaneous model
stress plus resolved covariance, including any remaining resolved variation.

Publish bounded sensitivities of incoming k/dissipation, rather than infer
an exact k from the unmeasured spanwise variance. For rough observational d,
use nominal 3.5 mm with declared 2.5/4.5 mm coordinate scenarios; do not move
the physical wall or change Yap distance. These are assumption ranges, not
experimental confidence intervals. Report source/digitization bounds,
numerical uncertainty pending CV4, parameter identifiability, common-fit
misfit and any conditional fits separately.

Freeze versioned S/R parameter manifests and native incident profiles before
reserved-flat checks. Preserve irreducible canonical-model misfit as a result;
agreement with every flat quantity is not required for milestone closure.
Do not mark a timed-out, undeveloped or untested candidate calibrated.

## Execution and retention

Use the private Linux source/build under
`/media/jay/新加卷/nektar-v2026-runs/cao-tamura-calibrated-validation`, and
store CV3 cases under `calibration/cv3/`. Leave the primary HIT checkout and
shared installation untouched. No more than 16 project CPU workers in total;
all OpenMP/BLAS threads are one. Every solve has an explicit allocation,
hard timeout, process/physical-time monitoring and command/exit receipt.
Large fields, histories and logs stay on the Linux volume. This protocol,
calibration manifests, compact data/evidence and figures are committed on
`codex/cao-tamura-calibrated-validation` when CV3's actual closure gate is met.

## Numerical preflight amendment (before parameter selection)

The geometric outer turbulence interpolation produced a smooth-floor epsilon
undershoot near the inlet at x=-0.958m,z=0.233m, with positive initial fields
and accepted weak-divergence residuals. The first valid field tending to zero
was retained; the rejection was not floored. Test streamwise/normal resolution,
finite-difference/linear solver controls, and the following explicitly inferred
outer-input scenarios before selecting the production discretization:

- Outer dissipation length 5mm versus the 25mm trial value; no outer
  dissipation was measured. This is an input diagnostic, not a source result.
- Arithmetic versus geometric C2 interpolation of the unmeasured k/epsilon
  between 0.8delta and 1.2delta. Both have positive endpoint inputs; neither
  clips or resets evolved fields. Record any effect on calibrated observables.

An unstable input/discretization is inadmissible for calibration. Choosing an
outer assumption for numerical admissibility must remain explicit and retain
its unresolved physical/numerical uncertainty until CV4; no hill targets enter
that choice. Original unsuccessful pilots remain part of the evidence.

A further input diagnostic applies a C2 upper saturation to the trial mean
log-polynomial before outer blending. The polynomial can slightly exceed
U_inf, introducing an artificial outer jet and large nondimensional CLS
strain where k/epsilon is large. The optional cap equals the polynomial up
to 0.98U_inf, equals U_inf from 1.02U_inf, and has matching first/second
derivatives at both joins. This constrains only prescribed inlet/initial
data, never evolved fields; fit residuals still use the unaltered frozen
observations. It is a declared inlet-family assumption to assess, not a
claim that the experiment measured an exact cap.

Numerical coupling diagnostics also compare the continuous Q3/Q1 and Q3/Q2
velocity/pressure pairs on the identical physical quadrature grid, and retain
a tight pressure-inverse tolerance while reducing only velocity/scalar native
preconditioner accuracy. These choices change the discretization/iteration,
not closure coefficients or wall continuity. Pressure degree and both current
operator/native inverse tolerances must appear in case/command receipts.

The native SBDF3/IMEX3 method is considered as a transport stability diagnostic.
Its explicit imaginary-axis amplification is compared with the current IMEX2
tableau before testing it on the precursor. This changes time discretization
only; validate manufactured temporal order, equilibrium recovery and full
history restart (including the SDIRK startup) before using it for calibration.
Preserve the default IMEX2 behavior. IMEX3 uses two inertia coefficients: the
constant-diagonal SDIRK startup and the subsequent SBDF3 stage.

An additional dissipation inference follows the measured/inferred inlet shear:
epsilon_in=sqrt((a_epsilon*k_in*|dU_in/deta|/Q_in)^2+epsilon_outer^2),
with the same declared H..N policy trace collar. It matches local production
and dissipation for the standard model at unit multiplier when the outer
term is negligible, and approximately the CLS neutral shear using the common
Q_in. The outer term ensures a positive prescribed free-stream endpoint; no
evolved field is floored. This tests the mismatch between a canonical-log
epsilon guess and the non-log measured inlet gradient. The analytic derivative
of the actual C2 mean inlet is used, not a numerical data derivative or hill
information. Treat this as a named inference family; retain the original
NeutralLog inference and all unsuccessful pilots.

IMEX3 did not resolve the coarse-grid precursor positivity rejection. Its
CLS velocity-order verification on the tested p6 mesh was also unaccepted.
Do not enable it in the delivered solver or use it to freeze CV3. The temporary
binary, native fields, source/binary hashes and failed checks remain diagnostic
evidence. IMEX1/2 remain the supported methods.

A discrete startup continuation is now tested to separate the abrupt
initial/inlet compatibility front from a developed unforced solution. Cache
the complete initial native weak residual R0, including wall/outlet loads,
and add w(t)R0 as a temporary load, with w the C2 fade from one to zero over
a declared duration T0. No evolved field is reset or clipped. Only an
experimental plane without external forcing may use this option. From
t>=T0 the additional load is identically zero and the original transported
models/wall equations apply. This warmup is not physical experimental time.
Require at least three unforced flow-through times after T0, then two small-
drift one-flow-time extensions. Record startup weight in every native stage.
Mid-warmup history restart is explicitly unsupported because R0 is not in
the history payload; post-warmup restart must pass the usual MPI/history test.
Before use verify residual cancellation at the initial field, the modified
weak residual, positive variables/continuity, exact zero forcing after T0,
post-warmup restart, and sensitivity to the warmup duration.

The cancellation test identified a generator bug: the prepared epsilon collar
used the selected policy's EpsilonWall, but the retained floor Dirichlet value
was still the neutral-log Ustar^3/(Kappa*H). VCS reapplied that value at startup.
With common physical K this produced a 0.0043467 (4.9%) epsilon jump in rough
original CLS before its first wall update. The floor Dirichlet initializer now
uses EpsilonWall for every policy, including varied K. Pre-fix pilots are not
admitted to the calibration design. The corrected hold control reduces the
first-map epsilon difference to approximately 3.2e-12.

The native continuity audit also accessed a boundary orientation-sign array
unconditionally. That array is unused when the pressure map has no sign
changes (including Q1). Follow ContField's conditional sign handling. Verify
open-boundary low-order-pressure audit agreement before accepting its results.

Startup continuation is a numerical initializer. Hill sessions set its duration
to zero; matched physical-parameter comparisons may differ in that initializer
only when both statistics windows are unforced. A continuous run containing the
specified post-warmup 3TF and two 1TF checkpoints is equivalent to separate
history-continuation segments; retain its checkpoint/history provenance.

A measured-rms inlet family is added before examining any objective scores.
The constant neutral k seed retains excessive turbulence across the outer
layer, while the allowed Fig.3c longitudinal rms decays throughout that layer.
Fit a weighted cubic to log(sigma_u/u_tau) versus eta/0.25m using positive,
resolved Fig.3c points only. Prescribe k=a_k*a_var*sigma_fit^2, with a_var=1
as a reference assumption on the two unmeasured variances, a constant k below
N, a positive geometric N..2N collar and the declared outer blend. Vary a_k
on the original bracket; do not call the template measured k or infer exact
spanwise variance. Reserved Fig.3d remains excluded.

For this family prescribe a convex blend of compatible neutral shear at H
and the mean-fit derivative at N, integrated analytically from N. The measured
mean fit is unchanged above N. This removes an abrupt near-wall mismatch
without changing any evolved equation or wall model. For dissipation use
`ShearLength`: sqrt((a_e*k*U'/Q_in)^2 +
(0.09^0.75*k^1.5/L_outer)^2). The reference 0.09 defines an inferred outer
length only and does not replace the CLS Cmu. Use a positive geometric H..N
epsilon collar that matches each policy's EpsilonWall. All physical inlet
values at and above N are common to the two closures. Independent native
expression evaluation, positivity and analytic-derivative checks are required.

Audit the full native momentum-stress transverse tangent for these prescribed
states, including nonlinear CLS terms, in both perturbation directions and
at two finite-difference sizes. Negative tangents are recorded as a possible
cause of difficulty, not silently repaired by changing closure coefficients
or used alone to claim an instability mechanism. Keep raw-input and revised-
input diagnostics. Numerical admissibility still requires accepted native
unforced evolution and the stated development gates.

For throughput/admissibility pilots also vary the number of prescribed outer-
blend mesh cells, outer normal layers, streamwise columns and dt. The revised
rms family removes the earlier constant-k outer cliff, so the ten blend cells
used for that diagnostic are not automatically required. Any selected coarse
calibration mesh is explicitly provisional: retain its strong divergence,
one-sided traces and field drift, and reopen the fit if CV4 changes incident
observables. Do not infer mesh convergence from a small weak residual.

A further numerical extension tests constant-coefficient grad-div stabilization:
add gamma*(div(U),div(v)) to the native momentum weak form, gamma>=0, with
its full vector coupling implicit, in the independent stage audit, and in the
pressure preconditioner. This is zero for the continuum incompressibility
constraint, adds nonnegative discrete energy dissipation gamma*||div(U)||^2,
and changes neither the turbulence closure nor integrated AWF strip inventory.
The normal reaction absorbs its wall-normal contribution; the matching
auxiliary tangential traction remains unchanged. Record gamma explicitly in
restart/configuration receipts. Default gamma=0 preserves existing paths.
Test gamma=0.01 m^2/s in the present analysis units, bracket its effect, and
verify native weak residuals, neutral-equilibrium invariance, zero recovery,
MPI/history identity and canonical transport smoke tests before admitting it
to calibration. Finite-mesh bias and strong-divergence convergence remain CV4
gates; gamma is not fitted to experimental or held-out targets.

The dt=0.05 explicit-transport pilots remain unaccepted even with grad-div.
Their native residuals do not justify accepting a failing timestep. Test an
optional fully implicit advection partition in addition to smaller explicit
timesteps. Keep native IMEX1/2 inertia/startup and explicit external forcing,
put -U.grad(phi) in the coupled nonlinear stage for every transported field,
and supply the same current convection to the AWF strip inertial predictor.
The AWF auxiliary forcing still uses its current physical convection. Moving
the bulk term without the strip contribution would be inconsistent and is
not admitted. Preserve explicit advection as the default, distinguish the
partition in restart signatures, verify standard/CLS manufactured order and
cross-rank full-history continuation, and retain every rejected pilot. This
changes the time partition, never the canonical constitutive/source equations.

The declared design will now be attempted using the measured-rms/ShearLength
family, common reference Q/kappa, implicit advection, gamma=0.01, Q3/Q2,
6 outer layers plus 3 blend cells, streamwise columns [-1,0,1.2]m and dt=0.05
in the declared analysis units. Use serial DirectMultiLevelStaticCond per case,
separate case CPUs and one controller, at most 16 project workers in aggregate.
Retain centre/axis failures rather than infer global model infeasibility from
a failed centre. First attempt higher-k, lower-ks and higher-epsilon probes
for an early admissibility check; selection still waits for all design attempts
to conclude. Each initial solve covers at least 3TF, followed by up to twenty
1TF extensions. A solve has a 7200s stage timeout. These are computational
budgets, not evidence of physical nonexistence. An exhausted/drifting/rejected
case remains ineligible. Native full weak residual, resolved and strip mass,
independent surface differentiation, and independently reconstructed BDF2
strip momentum identities are required in addition to field drift. Record
strong flux/divergence defects without claiming their convergence.

Both canonical manufactured cases recovered second-order implicit transport,
and native prefix-to-MPI continuations passed. Independent reconstructed flat
strip momentum balances agree with the native point residual to below 5e-16
in the initial implicit pilots; the pointwise residual itself is not zero on
the coarse mesh. A three-checkpoint smooth standard control at 8–10TF still
misses the strict epsilon drift criterion (0.203% then 0.145%, versus 0.1%);
it is not calibrated. Initial fixture-only regression rejections (serial
backend requested on MPI; a hold duration outside the generator's declared
range) were corrected and are not solver-physics failure evidence.

If the original design yields no developed CLS candidate, run smaller-step
controls before attributing failure to the constitutive law. In particular,
compare dt=0.01 with backward Euler (native IMEX1 plus implicit advection) for
the S/R centres. BDF2's extrapolated state can become negative during rapid
decay even when its saved states are positive. These are diagnostic controls,
not an automatic replacement calibration; first-order accuracy and inventory
checks are required before admitting a separate calibration recipe. Keep the
completed second-order design and its rejected candidates intact.

The first complete design produced 13 developed standard cases and no eligible
CLS case; no common parameters were selected or frozen. A solver-domain audit
then found that Newton rejected a signed intermediate *linear Picard image*,
even when its actual candidate k/epsilon remained positive. The residual
x-F(x) is well-defined for such an image: closure evaluations use x, and the
linear image is not an accepted physical field. Test the corrected domain:
require finite images, bound all physical Newton candidate updates, keep a
positive incoming first guess if its initial image is signed, and require
strict positivity at convergence and output. No state is clipped or floored.
Record nonpositive-image counts separately from physical guard activations.
The first design remains retained as the pre-fix attempt. Require a previously
rejected controlled case, full native residual/inventory checks and regression
verification before rerunning the same physical design with this correction.

Before interpreting CLS failures physically, test quadrature adequacy. Q3
with five GLL points integrates degree seven polynomials exactly; even a
constant-coefficient cubic gradient stress times the degree-two test gradient
can reach degree eight. The actual Cmu/k/epsilon factors are nonpolynomial.
Increase the common velocity/scalar/pressure quadrature grid to eight GLL
points, without changing approximation degree, mesh or closure coefficients.
Compare the same longest-surviving S/R cases and a rejected standard control;
this tests aliasing separately from h/p and physical-parameter changes. A
point-sampled negative constitutive tangent alone does not identify the cause
of an underintegrated solver failure. Retain the five-point campaign.

The rough eight-point control still rejected a nonlinear stage before 3TF.
Do not call quadrature increase a remedy. Before concluding the numerical
preflight, compare gamma=0.1 against 0.01 at fixed coarse discretization and
test a separate streamwise/normal refinement with four streamwise columns,
ten outer layers and eight-point quadrature. These controls distinguish
unresolved discrete divergence and spatial error from the retained
constitutive-tangent diagnostic; they do not establish mesh convergence or
authorize a coefficient change. The symbolic manufactured-case constructor
also requires serialization: its SymPy global assumption caches are not
thread-safe, although the resulting native solver executions may run in
parallel. A fixture-generation exception is not a physical solve rejection.

All of these controls rejected a nonlinear stage; none supplies a developed
CLS precursor. The independent parallel-shear stress-tangent oracle passes,
including the two physical sides of the max-invariant branch. It demonstrates
negative local transverse principal coefficients for some sampled positive
states even with retained molecular viscosity. This is a local constitutive
diagnostic, not proof that it caused every full-flow rejection or that no
canonical developed solution exists. Keep CV3 open: do not freeze a shared
calibration, unlock reserved targets, or replace canonical CLS to obtain one.
The separate [attempt record](cao-tamura-cv3-calibration-attempt.md) retains
the numerical failures and conditional standard-only design minima.

Final classification distinguishes the old rough-centre line-search rejection
from the old rough ks-high and smooth k-low linear-image positivity rejections.
Test those latter two CLS design points explicitly with the corrected Newton
domain, the unchanged design physics/numerics and the full development gates.
The centre alone does not exercise the old CLS image-guard stopping condition.
A successful targeted candidate would justify a new complete design recipe;
retain both outcomes and do not overwrite the first design.

Both targeted full-development cases rejected a stage: rough ks-high advanced
through t=3.15 (previously t=1.45), with signed linear images and positive
accepted physical fields; smooth k-low failed continuation after t=6.40.
Neither is developed. Add a short 35-step rough ks-high native regression
through the former stopping point, requiring signed-image branch coverage,
positive accepted scalars, strict weak/projected/mass gates and independent
strip inventories. Its successful short evolution is software verification,
never a calibration substitute. Keep the first complete design intact;
no successful developed CLS point justified repeating all 32 cases.
