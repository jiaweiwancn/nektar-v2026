# Calibrated validation against the Cao–Tamura steep-hill experiment

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: plan / computational experiment
- Origin Date: 2026-11-01
- Verification Status: CV1 source/data and CV2 implementation/readiness VERIFIED within recorded limits; CV3 IN PROGRESS under user-approved EP-transition-first sequence; R1 reference/standard controls pass but CLS composite-floor gates remain unmet; CV4–CV8 PLANNED
- Version Label: cao_tamura_calibrated_validation_v2_ep_transition_first

## Objective, branch and starting evidence

Determine how well the implemented standard k–epsilon and high-Re CLS URANS
models reproduce the smooth and rough steep-hill measurements of Cao and
Tamura (2006), after calibrating the undisturbed boundary layers independently
of the hill measurements. Quantify the effect of the length and placement of
the upstream transition from equilibrium-preserving AWF (EP-AWF) to original
AWF, together with numerical, calibration and sampling uncertainty.

All implementation and subsequent experimental validation for this plan belong
on **`codex/cao-tamura-calibrated-validation`**, starting at
`c48d9fb14afeca6caa80cf5b3f8ce0a34a02d2d7`. Commit and push after each completed
milestone. Initial planning commit `038bff7` recorded the design; the CV1
source/data closure below is separate from numerical or experimental validation.

The [terrain implementation record](hill-wall-geometry-transition.md),
[WG4 evidence](hill-wall-geometry-evidence/wg4.json),
[rank comparison](hill-wall-geometry-evidence/wg4-rank-check.json) and
[public terrain example](../../solvers/IncNavierStokesSolver/Examples/TerrainHill/README.md)
are the starting evidence. They establish geometric coupling, spatial masks,
short hill solves, restart/MPI checks and neutral reference controls. The hill
pilots cover only two timesteps, ending at 4e-5 s. Their nominal inlet,
roughness, roof and span are verification inputs, not experimental calibration.
The five periodic neutral controls do not replace inlet/outlet convergence
tests or long-time hill validation.

Questions to answer are:

1. Can separately calibrated smooth and rough approach flows reproduce the
   measured mean profiles, crest speed-up and downstream reattachment?
2. How do standard k–epsilon, high-Re CLS, compatible log law and the two AWF
   policies differ when inlet and surface calibration are held fixed?
3. Does the upstream EP-AWF transition leave a developed-flow effect larger
   than numerical and sampling uncertainty? Is there a width-insensitive range?
4. Which discrepancies arise from inflow, finite-Re transport, wall modelling,
   numerical error or missing three-dimensional dynamics?

Improvement by CLS or EP-AWF is a hypothesis, not an acceptance requirement.
A completed comparison can demonstrate that a model misses an experimental
target; that result must not trigger retrospective fitting to the hill.

## Experimental specification and separation of data

Primary source: S. Cao and T. Tamura, *Experimental study on roughness effects
on turbulent boundary layer flow over a two-dimensional steep hill*, Journal
of Wind Engineering and Industrial Aerodynamics 94 (2006), 1–19,
[doi:10.1016/j.jweia.2005.10.001](https://doi.org/10.1016/j.jweia.2005.10.001).
The local PDF is
`docs/literature/Experimental study on roughness effects on turbulent boundary layer flow over a two-dimensional steep hill.pdf`.
Its SHA256 at planning is
`adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df`.
Use the PDF as the source of figures, units and curve labels; the adjacent OCR
Markdown has unit/character errors. Do not add the literature collection to Git.

Use physical model scale throughout; denote hill height by **Hh**, and AWF
interface normal height by **hm**, so these two heights cannot be confused.

| Quantity | Experimental specification | Treatment |
|---|---|---|
| Hill | Hh = 0.040 m; half-length L = 0.100 m; yh(x) = Hh cos²(pi x/(2L)) for abs(x) < L, zero elsewhere | Same compact C1 shape in smooth and rough cases; crest x = 0, feet x = ±0.100 m |
| Facility | Test section 1.0 m wide, 0.8 m high, 7.0 m long; hill span 0.8 m with end plates; 5% blockage | Audit end-plate geometry and ceiling/sidewall conditions; a narrow periodic span is a proxy, not this facility |
| Fetch | Crest/reference station 5.0 m from contraction exit | First 1.0 m has decreasing-height turbulence-generating cubes in BOTH cases; smooth does not mean a laminar inlet |
| Surface | Remaining floor bare for smooth; 5 mm staggered cubes for rough, also on the rough hill | Frontal-area roughness density 0.041; geometric cube height is not equivalent sand-grain ks |
| Boundary-layer depth | delta about 0.250 m = 6.25 Hh for both inflows | Preserve finite depth and outer flow; do not extend a log profile to the roof |
| Friction velocity | u_tau/U_inf = 0.033 smooth, 0.051 rough | Flat-floor calibration constraints, not hill fitting parameters |
| Aerodynamic roughness | z0 = 4e-6 m smooth, 2e-4 m rough | Positive z0 even for the compatible smooth log wall; AWF uses separately identified ks |
| Displacement | Rough d = 0.7 times cube height = 0.0035 m | Assumed by the authors, not directly measured; retain this distinction in uncertainty analysis |
| Reynolds numbers | Re_H = U_H Hh/nu about 12000 smooth and 10500 rough; U_H is undisturbed velocity at hill height | Do not substitute U_inf in this definition; reconcile with rough Re_s = u_tau z0/nu about 3.97 |
| Inflow intensity | sigma_u/U_H = 0.0893 smooth, 0.151 rough at hill height | Additional flat-floor constraint/diagnostic, with the stress limitations below |
| Sampling | About 30 s at 1 kHz per measurement point | Correlated samples; not 30000 independent observations |

Recover U_inf/nu consistently from the source and digitized profiles before
choosing dimensional velocity and viscosity. If absolute temperature/speed is
unavailable, specify a nondimensional Reynolds-matched experiment and label
any dimensional reconstruction as inferred. Do not report the pilot's u_tau,
nu or 2 mm ks as measured experimental values.

**Calibration data:** Fig. 3(a–e), Table 1, flat-floor u_tau, delta, z0 and the
reported hill-height incident intensity. Preserve some heights or one profile
as flat-floor cross-check data. The paper does not provide measured epsilon
or the spanwise variance needed for a unique k profile. Document assumptions
and bounds instead of calling inferred k/epsilon experimental measurements.

**Held-out hill data:** mean profiles (Fig. 4), intensity (Fig. 5), fractional
speed-up (Fig. 6), crest turbulence (Fig. 7), wake mean/dividing points and
near-floor traverses (Figs. 8–9), wake rms (Fig. 10) and the available Reynolds
stress (Fig. 11). Extract actual station/height lists from the figures in CV1.
Reynolds stress inside the separation bubble was not measured with the
split-fibre probe; do not create such a validation target. Spectra and
reverse-flow intermittency are contextual diagnostics, not primary URANS
validation targets.

Primary endpoints are mean velocity, crest speed-up and reattachment.
Total rms and Reynolds stress are secondary closure diagnostics, compared
where the measurement and observation operator support them.

| Held-out scalar | Smooth inflow / smooth hill (SS) | Rough inflow / rough hill (RR) |
|---|---:|---:|
| Crest speed-up at vertical y' = 0.010 m = 0.25 Hh | 0.457 | 0.875 |
| Reattachment x_R/Hh, measured from crest | 5.4 | 6.5 |

The additional Fig. 6(b) crest targets are 0.918 for rough inflow / smooth hill
(RS) and 0.245 for smooth inflow / rough hill (SR), at the same height. Full
wake data are not assumed available for these mixed cases. Verify the physical
extent of rough/smooth coverage before specifying their material boundaries;
if the paper is insufficient, declare the ambiguity and bracket it.

## Models, walls and physical compatibility

Keep the canonical model coefficients and CLS physical closest-wall distance.
The primary route for BOTH SS and RR is EP-AWF upstream, a transition before
the hill, and original AWF downstream. Uniform-original-AWF and compatible-log
comparisons remain separately tracked. Successful all-original calibration is
not a prerequisite for primary EP calibration or its hill study. Retain the
ten requested model/wall/surface case slots with explicit validity/status:

| Closure | Compatible neutral log law | Original AWF | EP-AWF with fixed correction times spatial mask |
|---|---|---|---|
| Standard k–epsilon | Separate requested baseline | Separate requested baseline | Primary transition study |
| High-Re CLS | Not supported by this log policy | Separate requested baseline | Primary transition study |

Flow-switched EP-AWF is a secondary comparison for both closures and both
surfaces after fixed-mask behaviour is understood. Calibrate the primary
composite route first. Within a closure, transfer the same calibrated inlet
and physical surface parameters to original/log comparisons; record any
approach-flow mismatch rather than hide it by retuning. Each comparison needs
its own valid native evolution and observation evidence; it does not inherit
the primary route's calibration-complete status. Record failed/unexecuted
baselines explicitly while the validated primary EP study advances.
The log wall uses the corresponding calibrated positive z0, not AWF ks.
First report a common-surface-parameter comparison across closures. If distinct
closure-conditioned best fits are also useful, label them separately and
publish both; they cannot establish a pure closure effect.

Start the fixed correction at chi_fixed = 1 where admissible. Freeze and
report the failure policy, cutoff/high-Re checks and effective correction
coverage before production. `RequirePreservation` is appropriate for compatible
analytic verification. A measured-inflow case using `UseOriginal` must report
fallback fractions and locations; a case that falls back everywhere is an
original-AWF case, not evidence about active EP-AWF. Do not clip an inadmissible
cutoff or suppress the fallback receipt to complete a model matrix.

Use **retained resolved molecular viscosity** (`Molecular`) as the primary
finite-Re experimental track. Also run selected `HighRe` resolved-transport
variants, which omit that contribution, for sensitivity and the existing exact
neutral-reference tests. Both policies retain positive molecular viscosity
inside the AWF auxiliary problem. A high-Re turbulence closure can use a
smooth-wall AWF without resolving the viscous sublayer, provided hm is valid;
this does not introduce the low-Re CLS model.

Three compatibility issues must be settled before fitting or production:

1. **Smooth wall:** AWF ks = 0 is supported, but original-AWF smooth plane
   tests and rough-wall hill pilots do not validate smooth experimental hills.
   Verify the smooth EP-AWF matching anchor, admissibility and d = 0 reference
   on their own. Do not keep the 2 mm rough anchor while setting ks to zero.
2. **Rough displacement:** current AWF/CLS reject nonzero `RANSHeightShift`;
   standard log uses a nonnegative *additive* height shift. The paper's
   log(y-d) is not obtained by entering d in that parameter. Separate physical
   normal distance eta, data-fit displacement d and aerodynamic z0/ks. Do not
   move the geometric wall, change the CLS Yap distance to eta-d, or remove
   the physical 0..hm strip inventory to imitate displacement.
3. **Reference equilibrium:** derive/check the bulk and wall residuals for
   each proposed reference. Physical-distance CLS/Yap, finite nu, nonconstant
   k and finite-depth profiles need not admit a displaced constant-k analytic
   log equilibrium. The models' coefficient-derived kappa need not equal the
   measured log-fit slope. Exact d = 0 HighRe preservation remains a separate
   verification claim. Do not claim exact preservation of the measured finite
   boundary layer merely because the EP-AWF option is enabled.

CV2 must record the resulting decision. The default experimental route is a
measured/inferred finite-depth inlet and a computed flat-floor precursor whose
profile at the would-be hill location is calibrated to the data. Use the
current canonical equations where consistent, and expose remaining structural
misfit. If a displaced wall-reference extension is necessary, derive its
traction, dissipation, production, epsilon data, sensors, corner compatibility
and geometric inventory before implementation. Test d -> 0, admissibility,
zero-mask recovery and restart signatures. Any altered closure or Yap physics
is an explicitly named secondary model variant, with independent tests; it
must not silently replace canonical CLS. An incompatible analytic reference
is an honest negative derivation result, not a reason to tune model constants
against crest speed-up.

## Terrain coupling and observation operator

**Retain integrated strip continuity for curved-wall AWF.** With physical wall
arc length s, tangent t, into-fluid normal n, signed curvature c = dt/ds dot n,
J(eta) = 1-c eta and V(hm) = hm-c hm²/2, keep

```text
J(hm) W_hm = -hm d_s U_P - V(hm) d_span V_P .
```

Keep its P-height counterpart for sensors, native differentiation of pulled-back
fluxes, conservative normal multiplier/adjoint and vector implicit solve.
Only the compatible log-wall policy prescribes an impermeable interface.
The auxiliary AWF remains locally parallel; its omitted sublayer tangential
diffusion and normal-stress gradients remain modelling limitations.

Require 0 < P < hm < N, 0 <= ks < hm, unique normal sampling and the existing
curvature/Jacobian guards. For the paper hill, max abs(c) is about 19.74 1/m;
the default abs(c hm) <= 0.25 implies hm <= about 0.0127 m. Select hm, P and N
jointly with calibrated ks and the smooth-wall matching range. Do not relax
geometric guards to accommodate a fitted roughness. Report any lack of a
jointly admissible strip. The current 8/4/16 mm hm/P/N is only a starting trial.

The compact hill is C1, with curvature jumps at the feet. Existing AWF global
pointwise normal-error maxima do not converge there because a continuous
velocity trace faces different one-sided targets. Preserve the actual hill
shape and report strong L2 error, global maximum, both targets, their half-jump
bound and integrated conservation. Require L2/refinement and observable
convergence; weak residual alone is insufficient. Any smoothed-shape study is
a separate geometry sensitivity, never the main reproduction.

Build and independently check the experimental observation operator:

- Paper coordinates are streamwise x, vertical y, and **vertical** height
  y' = y-yh(x). Solver 3D coordinates are x/streamwise, y/span, z/vertical
  (solver 2D uses x,y). A point at fixed experimental x,y' is not a sample
  at normal height eta = y' on a sloping hill.
- Interpolate native fields at actual physical probe coordinates with
  deterministic MPI ownership. Rotate AWF-frame vectors/tensors to the
  laboratory axes. Compare U_x, vertical variance and x–vertical covariance;
  in 3D these involve R_zz and R_xz, not R_yy/R_xy. Check with manufactured
  polynomial/vector/tensor probes and rank changes; do not use nearest GLL rows.
- Reproduce S(x,y') = [U_hill(x,y')-U_floor(y')]/U_floor(y') with the
  matching-surface/inflow flat reference and the paper's reference-station
  convention. Also report an x-local flat reference as a drift diagnostic,
  without changing the experimental definition. Flag near-zero denominators.
- Measurements below hm need a declared AWF reconstruction operator; they are
  not resolved CFD values. Reconstruct auxiliary tangential velocity using the
  accepted instantaneous forcing/metric, transform frames and verify value/shear
  matching at hm and auxiliary P evaluation. Verify independent outer N sampling
  against the resolved field; N is outside the auxiliary strip. Reconstruct
  normal flow consistently with strip continuity if needed. Publish
  reconstruction uncertainty and distinguish it from resolved interpolation.
  Log-law extrapolation through separated flow is not a reliable substitute
  for a missing reconstruction. A mean-velocity reconstruction does not supply
  sub-interface model stresses automatically: derive/test a separate consistent
  stress operator or mark those rms/stress comparisons unavailable.
- The paper's near-ground traverses are y = 2 mm for smooth and 7 mm for rough
  (2 mm above the cube tops), x = 150..320 mm at 10 mm spacing. Both can lie
  below a typical AWF interface. Homogenized ks does not resolve local flow
  among individual cubes. Bracket the observation height and hm dependence;
  if the wall model cannot support the measurement, mark the comparison
  unavailable rather than present an invented resolved value.
- Reproduce the U = 0 dividing-point extrapolation used for x_R. Report the
  near-probe zero crossing and wall-shear zero as additional diagnostics;
  neither is automatically the paper's x_R. Distinguish x_R from the bubble
  length x_R-x_separation; the paper did not measure the separation location.

For turbulence statistics, use time-averaged instantaneous model stresses plus
the resolved URANS covariance, in the same laboratory frame:

```text
R_total_ij = mean(R_model_ij(t))
           + mean((U_i(t)-mean(U_i)) (U_j(t)-mean(U_j))).
sigma_i = sqrt(R_total_ii).
```

Do not add rms values, apply a closure only to the time-mean velocity instead
of averaging its instantaneous stresses, or include molecular viscous stress
as Reynolds stress. Record negative diagonal stresses/realisability failures;
do not clamp them to improve agreement. Standard two-equation isotropic stress
cannot generally match the measured longitudinal and vertical variances
simultaneously. With no measured spanwise variance, inferred k has an explicit
assumption range. URANS coherent spectra and sign changes do not establish the
full measured broadband spectrum or unresolved reverse-flow probability.

## Calibration, domain and case design

Calibrate separate **S** and **R** empty-floor boundary layers with the actual
EP-AWF -> upstream transition -> original-AWF treatment, at the location where
the hill would be placed. The matched floor retains the hill case's inlet,
surface, boundary conditions and exact correction mask. Use a finite-depth
profile with a smooth outer blend, positive k/epsilon and the measured low
outer turbulence, not zero
turbulence at the roof. Fit a small identifiable set: inlet profile/fetch,
rough AWF ks (smooth ks fixed zero), and bounded turbulence length/profile
parameters. Use the flat shear stress to constrain u_tau independently of
mean-profile fitting. Report inferred epsilon and any spanwise-variance
assumption; use a range rather than an unsupported exact k measurement.
Keep canonical C_mu/CLS/Yap constants and AWF smooth constants fixed.

Use profile-height weights or fixed height bands so digitization density does
not set the fitting objective. Report weighted mean-profile error, stress
error, delta and incident intensities separately, parameter sensitivities,
identifiability and residuals. Initial engineering goals are U-profile error
<=3% of U_inf and u_tau/delta errors <=5%, subject to source precision; these
are planning targets, not established experimental confidence intervals.
If the canonical closure cannot meet all flat-flow quantities, preserve the
misfit and propagate it to hill interpretation. Do not compensate with a hill
fit. Freeze a versioned, policy-specific primary calibration manifest before
inspecting validation scores; later justified calibration revisions get a new version and rerun
affected cases, with earlier scores retained.

The user-approved [CV3 protocol v2](cao-tamura-cv3-calibration-protocol.md)
supersedes the [original-first v1 protocol](cao-tamura-cv3-calibration-protocol-v1.md).
The immediate task is **CV3-R1: reference recovery and composite-floor
feasibility**. Repeat the accepted CLS EP-AWF reference on the current solver,
then establish four valid unfitted composite precursors (S/R x standard/CLS)
before a new parameter search. Use the provisional width 2Hh=0.080 m and
endpoint xb=-0.200 m, so xa=-0.280 m; record this mask in the primary calibration
anchor. CV5 still determines width/placement sensitivity. Freeze the validated
primary policies and mark baseline transfers separately; uniform-original
success is not a freeze prerequisite. The v1 original-first campaign/selector/
freezer remain legacy tools until adapted/versioned for this sequence.

Start with a central-plane 2D URANS campaign. Match the finite-depth inflow
and 0.8 m facility ceiling in the primary geometry; assess roof boundary
conditions and sidewall/end-plate blockage from the source. Truncated inlet
locations x/Hh = -15 to -25 and outlets x/Hh >=20 to 30 are starting domain
candidates, not accepted sizes. Provide sufficient full-EP approach length
upstream of the widest mask. Check inlet/outlet placement independently; the
last measured profile at x/Hh = 7.5 is not a justified outflow location.

Use selected fully 3D spanwise-extruded cases to test 2D assumptions. Periodic
spans, initially considering L_span/Hh = 2, 4, 8, need span/resolution studies
and remain proxies for the finite end-plated tunnel. Match or bound the actual
end-plate/sidewall effect before claiming facility equivalence. General
non-extruded terrain and explicit resolution of the cube array are outside
this plan. Mark an unresolved facility effect as a validation limitation.

For EP-AWF, use the delivered C2 envelope at physical wall-foot x:

```text
r = clip((x-xa)/(xb-xa), 0, 1)
w = 1 - 10 r^3 + 15 r^4 - 6 r^5
chi = w * chi_fixed             or             chi = w * chi_flow .
```

The same chi must select strip dissipation and epsilon at P/hm; original
traction and production remain. Require exact original-AWF recovery at zero
mask. Use Plane for flat controls and CosineHill for hills; legacy Flat does
not support the spatial configuration.

| Width Delta/Hh | Delta (m) | xa (m), with xb = -0.200 m |
|---|---:|---:|
| 1 | 0.040 | -0.240 |
| 2 | 0.080 | -0.280 |
| 4 | 0.160 | -0.360 |
| 8 | 0.320 | -0.520 |

The downstream endpoint is strictly upstream of the foot at -0.100 m. The
mask is zero over the hill and wake. Keep xb and the mesh identical while
varying width, and resolve the narrowest ramp. Study placement separately
at a fixed width, initially xb = -0.240, -0.200, -0.160 m, all upstream.
Material roughness changes in RS/SR are independent of this correction mask;
align material jumps with wall elements and retain one-sided property traces.

The fixed-mask width matrix is two closures × two inflows/surfaces × four
widths × hill/empty-floor = **32 developed primary cases**, with original-AWF
and standard-log comparisons on separate tracks. Failed/unexecuted baselines
do not certify those comparisons or prevent valid primary cases advancing.
Cases reused by the final validation matrix keep one
case ID and fingerprint; do not rerun or count them as independent evidence.
Use flow-switched, molecular-policy, displacement, placement and 3D variants
as explicitly selected sensitivities rather than a full uncontrolled Cartesian
product. Publish the selection before examining hill errors.

For every width/placement, record floor-only drift and incident U, k, epsilon,
stress, boundary-layer thickness, chi and all sensor/fallback states at common
stations. Compare hill-minus-matched-floor changes as well as raw hill fields.
Select a usable width from approach-flow admissibility, numerical resolution
and a plateau, not from which width fits crest speed-up best. A provisional
plateau requires two adjacent widths to change incident mean velocity by
<1% U_inf and incident stress by <2% u_tau², with differences also assessed
against uncertainty. Numerical/sampling uncertainty in the differences must
be small enough to resolve those bands; an effect hidden by larger uncertainty
does not establish a plateau. Report S, x_R and stress sensitivity for all
widths even if no plateau exists; then widen the range or report unresolved
dependence.

## Numerical verification, stationarity and uncertainty

Carry forward the existing kernel/solver acceptance gates: relative kernel
errors <=5e-10; stage weak-divergence/projected-normal checks <=1e-7 in the
verification cases; full-field MPI/history differences <=1e-6 scaled by
max(1, field magnitude). Check primary and derived fields, positive k/epsilon,
true nonlinear/linear residuals, inactive-correction recovery, integrated
normal flux and mass/source budgets. Record production residuals in both raw
units and stated U_inf/Hh scaling; do not silently reinterpret a dimensioned
verification tolerance. Guard activation is a diagnostic to resolve, not a
successful case hidden by field clipping.

Numerical studies must include:

- Three successively refined **whole-domain** meshes at fixed hm/P/N, refining
  the ramp, both feet, separated shear layer, outer layer and wake. A first-cell
  change alone is not a spatial convergence test. Keep geometry interpolation
  and partition precision sufficient on every level.
- Separate polynomial-degree and three-timestep sequences (dt, dt/2, dt/4),
  with adequate quadrature and the same physical sampling windows. Separate
  admissible hm/P/N studies from mesh refinement; do not conflate a wall-model
  height change with discretization error.
- Domain/roof/outlet sensitivities, selected span refinement, and tolerance
  studies for representative narrow/wide-mask cases. Neutral preservation
  errors need their own inlet/outlet h/p/dt refinement, reporting U/k/epsilon,
  drift rates and L2/Linf norms; periodic stationarity alone is insufficient.
- Quantify changes in profiles, S and x_R, not only residuals. Use Richardson
  extrapolation/GCI only in an evidenced asymptotic regime with appropriate
  refinement ratios. For nonmonotone data give an uncertainty envelope and
  extend refinement if a scientific comparison remains unresolved.

Initial numerical-uncertainty goals are <=1% U_inf for mean velocity, <=0.02
absolute S, <=0.1 Hh for x_R and <=5% of the relevant rms scale. These are
decision targets to check against attainable resolution, not measured results.
Refine representative smooth/rough, standard/CLS and narrow/wide cases, then
document why their error bounds transfer to other cases. Expand testing when
that transfer fails. Do not infer all-case convergence from one easy flat case.

Before the large campaign, profile the terrain solver on <=16 project CPUs.
The WG4 8-rank CLS two-step run took 3296.49 s under shared-server load;
452 native inverse calls hit their iteration cap at a requested 1e-14,
although physical/true-residual and MPI gates passed. This is not a throughput
prediction for a developed experiment. Separate native preconditioner accuracy
from the physical residual target and improve/tune the former with controls;
do not weaken physical gates or smooth the physical CLS max(S,Omega) branch
to gain speed. Benchmark time per accepted developed-flow step, memory and
checkpoint cost before budgeting CPU-hours. The earlier growing long-time
3D CLS mode is not resolved by a two-step hill run: include perturbed long-time
flat/hill controls and timestep/rank checks before interpreting 3D statistics.

For each accepted case use a fixed timestep and record physical time from
startup, including pressure adjustment. As an initial sampling allocation,
discard at least 10 flow-through times T_F = domain_length/U_inf, then collect
at least 20 T_F; extend until stationarity and uncertainty criteria hold.
For unsteady flow also cover at least 20 dominant cycles where identifiable.
Track flux, drag, x_R and representative wake probes; compare three consecutive
windows for secular trends. Use at least ten blocks with block length exceeding
twice the estimated integral correlation time, extending the run as needed.
Compute autocorrelation-aware/block-bootstrap confidence intervals; 30 s of
the experimental record does not prescribe an identical CFD sample count.
Initial sampling goals are 95% interval half-width <=0.5% U_inf for mean U,
<=0.1 Hh for x_R and <=5% of rms. A stationary steady URANS solution needs no
invented cycles; its modelled turbulence must still enter the stress comparison.

Publish signed bias, profile errors and station-by-station plots, using fixed
normalizations/ranges from CV1; do not pool away local reverse-flow failures.
Define profile RMSE = sqrt(sum(w_i (q_CFD_i-q_exp_i)^2)/sum(w_i)), with fixed
height-band weights. Provisional engineering accuracy targets are mean-profile
RMSE/U_inf <=0.05, crest S absolute error <=0.05, x_R error <=0.5 Hh and
rms-profile RMSE/q_rms_ref <=0.15, where q_rms_ref is the weighted rms of the
experimental sigma profile on that station's frozen comparison range. Report
the rms error scaled by U_inf as well; do not use that larger scale to hide
relative turbulence error. Freeze or justify revised scales before viewing
validation scores. These targets do not substitute for experimental uncertainty
or force a model to pass.

Keep numerical, calibration, digitization, sampling and facility/model-form
uncertainty separate. Do not combine bounds by root-sum-square without an
independence argument. The reported split-fibre pitch calibration deviation
is not a complete 95% measurement uncertainty. Give observable-level outcomes:
supported within stated uncertainty, disagreement, or inconclusive because
uncertainty/data/observation support is insufficient. Numerical milestone
closure requires trustworthy evidence, not universal agreement with experiment.

## Milestones and closure evidence

CV1–CV2 are **CLOSED**; CV3 is **IN PROGRESS**; CV4–CV8 remain **PLANNED**. Maintain closure notes and compact
receipts under `docs/plans/cao-tamura-validation-evidence/` as they are produced;
store large simulation outputs only in the Linux run root. CV1's data/evidence
paths now exist; later milestone outputs remain planned. Each closure records
tested revision/source hashes, commands,
case manifest, configuration, failures, outcome and limitations, then a commit
and push on this branch. Reopen affected downstream assessments after any
physics, calibration, observation-operator or numerical change.

| Milestone | Dependencies | Deliverables and closure gate | Status |
|---|---|---|---|
| CV1 — Source and data freeze | Planning commit | PDF/figure/page audit, geometry/units/station tables, independently checked digitized curves and uncertainty/provenance, calibration/holdout split, dimensional/Re consistency, unresolved facility/mixed-surface details explicitly recorded. Freeze versioned data and score definitions before fitting. | CLOSED |
| CV2 — Wall/reference and probe readiness | CV1 | Parameterized experimental generator; smooth-anchor checks; displaced-reference derivation/decision; any necessary declared implementation; native physical-point interpolation and below-hm reconstruction with frame/rank/continuity tests. Canonical d=0, finite-nu, original recovery, curved continuity, history and legacy checks pass on Linux. Record unsupported observables rather than manufacture them. | CLOSED |
| CV3 — Separate S/R calibration | CV1–CV2 | R1: current-solver EP reference recovery and four developed primary composite-floor cases (S/R x standard/CLS). R2: fixed-coefficient common/conditional fits, uncertainty and incident residuals, with four validated policy-specific S/R manifests. R3: reserved flat crosschecks after freeze. Original/log transfers have separate statuses; all-original success is not a primary freeze prerequisite. No hill scores enter selection. | IN PROGRESS; [R1 execution](cao-tamura-cv3-r1-reference-feasibility.md) |
| CV4 — Numerical and long-time readiness | CV2–CV3 | Profiled solver and bounded runtime estimate; three-level h and dt evidence plus p/hm/domain sensitivities; inlet/outlet preservation convergence; perturbed 3D stability checks; demonstrated stationarity/statistical procedure, mass/strip budgets and MPI/history equivalence. No large production matrix until throughput and observable error are defensible. | PLANNED |
| CV5 — Developed transition study | Primary CV3 freeze and CV4 readiness | Fixed-mask 32-case primary hill/floor matrix, upstream-only width and separate placement study, incident diagnostics, switched subset and selected resolved-viscosity sensitivities. Original/log baselines proceed separately where valid. Quantify numerical/sampling significance and select/freeze width by approach-flow/plateau rules, retaining every width and any non-plateau result. | PLANNED |
| CV6 — Locked SS/RR validation | CV1–CV5 for primary EP cases; each baseline's own readiness | Four primary standard/CLS x S/R EP-transition hill cases, reusing identical CV5 cases; original/log and switched comparisons where valid. Report all ten requested case slots with successful/failed/unexecuted statuses. Held-out profiles/S/x_R/stress require each case's verified observation operator, calibration provenance and uncertainty. No hill-driven retuning or inherited baseline-validation claim. | PLANNED |
| CV7 — Mixed surfaces and 3D checks | CV4–CV6 | RS/SR crest comparisons with audited material coverage, selected span/domain/3D checks for S/R and both closures, total-stress statistics and limitations. Bound or explicitly leave unresolved end-plate, homogenization and 2D effects; no unsupported mixed-wake or broadband/intermittency claim. | PLANNED |
| CV8 — Reproducible assessment | CV1–CV7 | Final observable-level accuracy/uncertainty tables, common-calibration versus conditional-fit interpretation, width/placement recommendation or unresolved range, convergence figures, case/data/build hashes, executable reproduction commands and restart instructions. All outcomes and failed attempts retained; working tree reviewed, milestone committed and pushed. | PLANNED |

The [CV3 v1 attempt record](cao-tamura-cv3-calibration-attempt.md) reports the
completed original-first 32-case bounded design: 13 developed standard cases
and no eligible CLS case, including the later numerical controls. Those results
remain intact. They do not execute or reject the revised composite-first route.
Primary S/R calibration manifests and reserved-flat crosschecks remain
unfinished. [CV3-R1 execution](cao-tamura-cv3-r1-reference-feasibility.md) records
the v2 reference, implementation and composite-floor feasibility checks;
those checks do not close the remaining R2/R3 calibration and crosscheck work.

CV2 may conclude that exact displaced equilibrium is unavailable for a
canonical model; CV3 may expose irreducible flat-data misfit. Those findings
must narrow subsequent claims or trigger a documented plan revision. They do
not authorize a silent new closure. If admissibility or numerical uncertainty
prevents a specified production comparison, retain it as unexecuted with the
reason and revise the scope explicitly rather than mark it validated.

## Implementation layout and reproducibility

Extend the existing workflow rather than edit nominal pilot values in place:

- `solvers/IncNavierStokesSolver/Examples/TerrainHill/`: add a documented,
  parameterized calibrated-experiment entry point and example manifests.
- `solvers/IncNavierStokesSolver/Tests/RANS/`: add calibration, physical-probe,
  campaign, statistics and comparison utilities; reuse geometry/initialization,
  full-field audits, restart and evidence infrastructure where applicable.
  The existing `generate_terrain.py` and `terrain_initial.py` are verification
  tools, not already calibrated experiment executables.
- `library/SolverUtils/Turbulence/` and the terrain/velocity equation-system
  files: change only if CV2/CV4 establishes a specific physics or solver need;
  include regression and public option/restart documentation with that change.
- `docs/plans/cao-tamura-validation-data/` (planned): small derived tables with
  source page/figure, series, axis transforms, units, extraction version and
  uncertainty. Retain raw digitization locally; publish derived material only
  with appropriate source attribution and redistribution rights.
- `docs/plans/cao-tamura-validation-evidence/` and
  `docs/plans/cao-tamura-validation-figures/` (planned): compact JSON/Markdown
  receipts and Linux-generated scientific PNG/PDF figures. No checkpoints,
  large meshes or full solver logs in Git.

CV1 records all dependencies. CV2 must deliver runnable prepare/calibrate/run/
compare commands with explicit flags and versioned manifests; they do not
exist merely because they are specified here. Use existing external RANS
Release-build support as documented by the terrain example. Bind source and
binary hashes, full-precision SpatialDomains library, Nektar/MPI/BLAS versions,
mesh/geometry/quadrature, inlet data, surface parameters, physical coefficients,
viscosity policy, hm/P/N, mask/sensors, timesteps, MPI/thread allocation and
observation/statistical definitions to each case ID. Hash the digitized data
and calibration manifest. Save actual commands, stdout/stderr and exit receipts.

Verify that new d/reference/property inputs enter history signatures. Keep
complete native solution history, physical time and the same-case restart
checks across ranks. A changed timestep/configuration must not reuse incompatible
history; start an explicitly documented new segment and establish stationarity
again. Never reset startup time and label a transient as developed data.

## Linux execution and resource allocation

Source/Git work uses `/Users/jay/Documents/nektar-v2026`. Linux is authoritative
for builds, solver tests and numerical/performance conclusions. At execution,
inspect revisions and working trees on both machines before synchronization;
preserve other work. Use an isolated source/build/install under
`/media/jay/新加卷/nektar-v2026-runs/cao-tamura-calibrated-validation/`, with
`tmp/`, `data/`, `calibration/`, `verification/`, `cases/`, `statistics/` and
`failed/` subdirectories. Keep all heavy files there. Reuse the installed
Nektar environment as a dependency, while ensuring the intended private solver
and geometry library are actually loaded. Leave the primary Linux checkout,
unrelated HIT runs and shared installation unchanged.

Connect by `ssh linux-usb`, falling back to `ssh jay@100.105.247.110`.
This project's earlier explicit limit remains **16 CPU workers total**, within
the repository's server-wide 72-worker limit. Set `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and `BLIS_NUM_THREADS=1`.
Count concurrent builds, MPI ranks, test workers and postprocessing together;
for example one 8-rank solve plus one 4-rank control leaves at most four
workers for everything else. Choose allocations from current server load,
including physical-core/SMT overlap with other jobs; do not reuse old affinity
lists blindly. Set compiler/test/MPI parallelism explicitly, never from all
reported CPUs.

Only after CV4 benchmarks, publish a CPU-hour/disk budget by case family and
use checkpointed, bounded-walltime run segments. Monitor process liveness,
physical timestep progress, iteration counts, residuals, guard events, disk and
checkpoint freshness; give each command a recorded hard timeout and retain
failed/incomplete attempts. Resume from compatible history with a new receipt,
not an unrecorded retry. A throughput failure is a reason to diagnose or reduce
the declared campaign, never to count a short transient as experimental
validation. This plan starts no builds or experimental runs.

## CV1 closure (2026-11-01)

CV1 is closed as a **source audit and data freeze**. No CFD case was fitted or
run, and no experimental-accuracy claim is made. Implementation/reference,
probe reconstruction, calibration and numerical-readiness work remain CV2–CV4.
Integrated strip continuity and the 16-worker project limit remain unchanged.

Delivered [frozen data v1.1](cao-tamura-validation-data/README.md), its
[snapshot hashes](cao-tamura-validation-data/snapshot.json),
[figure inventory](cao-tamura-validation-data/figure_inventory.csv),
[score definitions](cao-tamura-validation-data/scores.json) and
[CV1 receipt](cao-tamura-validation-evidence/cv1.json). The source PDF's
19 pages and its identity/geometry/units were audited against rendered pages;
publisher metadata confirms the DOI/title. The fourteen-figure inventory
distinguishes calibration, validation, duplicate and context-only material.
The dataset has 39 reported scalar facts and 1,151 digitized points across
98 plotted series. Eighteen legend-occluded symbols are retained with zero
weight; 1,133 points are marked visible. Repeated drawings and repeated plots
do not count as independent measurements.

Final review corrected the unpublished v1 draft's classification of four
Fig. 6a upwind/downwind speed-up series: they are additional validation data;
only its crest series duplicate Fig. 6b. The v1.1 revision retains all numerical
values and bounds, records the role change and preserves the draft locally.

Independent pdfplumber and pypdf parsing agrees on all selected symbol centres,
shapes and fills. The separate Poppler raster check passes on 472 isolated
symbols, maximum centre difference 0.179 PDF points against the fixed 0.5-point
extraction bound. All seven plot-page overlays were visually inspected. This
is independence of software backends, not of human readers or experiments.
Fresh extraction reproduces all three numerical extraction products byte for
byte. The audit checks coordinates, fixed height weights, source/recipe/snapshot
hashes, calibration/holdout separation and reported scalar identifiers. Four
negative checks reject a changed source PDF, altered frozen table, hill data
introduced into calibration (even with an updated file hash), and overwrite of
the sealed snapshot. Commands, versions, counts and hashes are in the receipt.

The data audit changes several assumptions for subsequent work:

- Figures 4/5/8/10 use global y/Hh; Figure 3c/d uses (y-d)/delta and Figure 7
  uses yprime/delta. Fig. 11's numerical delta75 is absent, so its dimensional
  height is not invented. The CSV/series metadata preserve these distinctions.
- Fig. 3a gives U_H/Uinf about 0.79351 (S) and 0.68670 (R). Together with the
  reported Re_H these imply Uinf/nu about 378067 and 382261 per metre; the
  resulting rough Re_s=3.899 differs from the reported 3.97 by about 1.8%.
  Absolute Uinf and nu are not separately specified. Freeze Reynolds-matched
  dimensionless comparisons; dimensional Fig. 9/12 amplitudes remain conditional
  on independently justified scale, with no hill-data scale fitting.
- Fig. 3e, converted with the reported friction velocities, differs from the
  repeated Fig. 3a profile by up to 0.05014/0.04351 Uinf (S/R). This unresolved
  source normalization/plot discrepancy is explicit. Fit Fig. 3a/b/c, reserve
  Fig. 3d, and keep Fig. 3e as a duplicate consistency check.
- The supplied Fig. 8 labels `2.53` and `3.0` are qualified as inferred 2.5
  and 5.0 stations using geometry and repeated profiles. Text/figure traverse
  extents differ. Reported ground-extrapolated x_R is distinct from finite-
  probe zero crossings. These discrepancies and the unspecified extrapolation
  method enter observation/uncertainty decisions rather than silent corrections.
- End-plate dimensions, roof/sidewall boundary-layer data, exact mixed-case
  downstream roughness coverage and a full experimental uncertainty budget
  remain unspecified. The recorded 5% blockage versus 4% full-section frontal
  area ratio also needs a facility interpretation. CV4/CV7 must bound or declare
  these limits. Spanwise variance and epsilon are not measured inlet data.

CV1 checks ran as serial, lightweight source/data audits on macOS 14.6 arm64,
Python 3.12.14 (plus the system-Python standard-library check), pdfplumber
0.11.9, pypdf 6.10.0, Pillow 12.3.0 and Poppler 26.05.0. Tested base is
`038bff7` plus the CV1 files identified by their receipt hashes. No compiler,
MPI rank or Linux solver was used. Both Linux SSH endpoints timed out during
readiness checks; Linux availability is not a gate for this paper/data-only
milestone. No Linux build, solver, stability or performance result is inferred.
Local PDF renders and development-attempt diagnostics are outside the checkout
under `/tmp/nektar-cv1`; the source-bound auditor reproduces them. The Git
delivery contains derived numerical tables and compact evidence, not the PDF
or its page images.

## CV2 closure (2026-11-01)

CV2 is **CLOSED** as implementation and numerical/observation readiness.
The [reference and observation decision](cao-tamura-cv2-reference-observations.md)
contains the displaced-reference derivation, inlet assumptions, operator
definitions and limits. The [Linux receipt](cao-tamura-validation-evidence/cv2.json)
binds source/report/binary hashes and commands. CV1's sealed data are unchanged.
No developed precursor was calibrated, no experimental hill score was inspected,
and stationarity, convergence and experimental accuracy remain CV3–CV8 work.

The canonical models and physical Yap distance are retained. A displaced
constant-k profile has unforced momentum residual
`-nu*u_tau/(kappa*(eta-d)^2)`; its epsilon residual additionally contains
`2*nu*A/(eta-d)^3` and the CLS physical-distance Yap difference. Thus a pure
shifted finite-nu log law is not an exact reference, and displaced CLS generally
also fails equilibrium with zero resolved nu. There is no new displaced closure
variant. The experimental route uses a finite-depth inlet and computed floor
precursor, with d observational and `RANSHeightShift=0`.

Delivered implementation:

- `Examples/TerrainHill/cao_tamura.py` provides generate/prepare/probes/run/
  observe/statistics/calibrate/compare commands and provisional S/R JSONs.
  Configuration covers model/wall, molecular policy, H/P/N, ks, fetch/domain,
  mesh/order/span, timestep, outer turbulence assumptions and upstream mask.
  Identical meshes are retained across the four mask widths. The trial rough
  ks=6 mm is unfitted. The nominal inlet U_H=1 dimensionalization fixes input
  Re_H; CV3 must verify the developed incident U(Hh) and actual Re_H as well.
- Native initial preparation evaluates the actual named inlet and its
  streamfunction primitive, checks positive Jacobians and independently audits
  velocity/scalars, positivity and span invariance. Positive geometric outer
  turbulence blends and adequate outer mesh layers avoid projection undershoot.
- `RANSObservationAudit` interpolates quad/full-hex native fields at laboratory
  physical points, including deterministic MPI traces and curved boundaries.
  Read-only AWF reconstruction retains the original auxiliary momentum profile
  and integrated strip-continuity normal flow. U_P and U_A(P) are distinguished;
  H value mismatch is reported. Nonfinite native inverse seeds/residuals are
  handled explicitly, fixing a failure exposed by the new floor geometry.
- Physical reference probes, speedup, covariance ratios and two-/three-tail
  ground-extrapolated x_R operators are available. Independent synthetic
  fixtures test these methods and reject incorrect vertical reference heights.
  Comparison requires a completed, versioned CV3 calibration before accessing
  reserved or hill scores. Flat-only candidate objectives are runnable now.
- Statistics use physical-time weights and mean model stress plus resolved
  covariance. Native configuration/probe/binary receipts bind their inputs,
  including compatible history segments. Prepared/zero-time/manufactured and
  unbound data cannot enter experimental statistics. Actual experimental inlet
  expressions/parameters and observational d enter history compatibility;
  legacy sessions retain their previous signature.

Authoritative Linux checks used base `ba3e005cde9d6d7dd7ab0f86e0ab31f0b69d4482`
plus the source contents identified by the receipt. The isolated source/build
and every heavy output are under
`/media/jay/新加卷/nektar-v2026-runs/cao-tamura-calibrated-validation`.
Tailscale SSH `jay@100.105.247.110` was used. The primary Linux HIT checkout,
shared build and installation were preserved. Compilation used at most 12
workers, native tests 1/2/4 MPI ranks, and all OpenMP/BLAS threads were one;
the combined project allocation never exceeded 16 workers.

| Gate | Linux result |
|---|---|
| Registered regression tests | All 11 CTests pass; the final added shear check also passes its targeted CTest. |
| Reference/reconstruction kernels | 975 checks, including 32 d/nu/model/height residual cases; auxiliary H shear-flux error <=8.52e-13. Finite-difference bulk residual checks use their declared 2e-8 scaled tolerance. |
| Smooth anchors | Both closures use ks=0 anchors recomputed from J0; largest three-step normalized reference error <6.8e-9. |
| Canonical neutral controls | Five models/wall modes, five flow clocks each: U error <=4.51e-8 and epsilon error <=6.74e-7; correction stays on its reference plateau. |
| Native points/frame/MPI | 2D/3D, horizontal/inclined polynomial/vector/tensor probes: error <1.93e-12; independent auxiliary quadrature <1.91e-12; curved rank difference <1.43e-14. Outside and unsupported probes are explicit. |
| Generator/preparer | 19 case definitions, 15 positive native preparations including a full-hex rough hill, identical width meshes and invalid-input checks. Full-hex input audit <2.92e-13, scaled span variation <3.73e-11. |
| Finite-nu and recovery | Both closures' Molecular controls converge with positive fields; zero-weight EP recovers original AWF within 1e-6. Exact finite-nu pure-log preservation is not claimed. |
| Curved continuity | Nonzero transpiration about 3.65e-3 in the verification units; projected normal <=1.80e-10, weak divergence <=5.83e-10, actual normal-plus-side flux defect <9.68e-14. |
| Experimental workflow/history | S/standard-EP and R/CLS-EP two-step Molecular floor solves; no guard activations, weak divergence <1.10e-10, mass defect <1.77e-11. Two-to-one-rank history differences <8.78e-7 over primary/derived fields; eight changed-input continuations rejected. Cross-history statistics pass. |
| Legacy acceptance | 20 standard-AWF cases, independent wall/kernel checks and three high-order CLS pilot cases pass, including original recovery, budgets and 1-to-4/4-to-2 rank history. |
| Observation/statistics oracles | Time/covariance error <1.74e-18; negative diagonals/missing coverage retained; speedup and ground extrapolation exact on independent fixtures; wrong height and pre-freeze score access rejected. |

The C1-foot limitation remains. In these short curved verification cases the
largest target jump is 2.771e-3, with half-jump bound 1.386e-3 and measured
maximum normal error 1.787e-3. These are raw verification units, not calibrated
experimental errors. Earlier strong-L2/refinement evidence is retained; this
milestone does not claim new maximum-norm convergence at the feet. Exact-foot
sub-interface probes are unavailable. Auxiliary momentum-gradient continuity
is also not imposed: the checked ODE shear identity and HighRe reference
traction must not be confused with exact gradient matching of the outer model.

Below-hm k/epsilon/full Reynolds tensors remain unavailable, as do log-wall
sub-interface values and auxiliary reconstruction for unsupported body forcing
or dealiasing. Current flat stress/rms objectives therefore report incomplete
height coverage rather than a complete-range score. Unknown experimental
delta75, unscaled dimensional amplitudes, separately averaged wall-shear zeros,
facility effects and actual tail/grid/sampling uncertainty retain their stated
conditions. The new seed fits only allowed Fig. 3a data; inferred turbulence
and the source displacement assumptions are not new measurements.

Development failures are retained in the run root: initial compiler/test-driver
errors; curved-segment inverse initialization; coarse outer projection
undershoot; incorrectly supplying observational z0 as an AWF policy parameter;
nonfinite terrain inverse seeds; and stale expected test counts/rejection text.
Final suites use `controls-e`, `probes-f`, `generator-c`, `workflow-e`,
`utilities-d`, `preservation-final` and `legacy-final-b`, with full logs/fields
outside Git. The unchanged solver binaries were reused while final independent
postprocessing and the targeted interface-shear kernel check were added.
