# CV2 reference, inlet and observation decisions

This records the implementation decisions for CV2 of the
[calibrated-validation plan](cao-tamura-calibrated-validation.md). It does not
constitute a developed flat-floor calibration or experimental validation.
The frozen inputs remain `cao-tamura-2006-v1.1` from CV1.

## Canonical displaced-reference residual

Let eta denote **physical normal wall distance**, d the observational
displacement, r = eta-d > 0, and define

```
K = Q u_tau^2,   A = u_tau^3/kappa,
U = U_H + u_tau/kappa log(r/r_H),
epsilon = A/r,   nu_t = kappa u_tau r.
```

Q and kappa are the implemented coefficient-derived neutral references, not
the experimental log-fit constants. For standard k–epsilon,
Q = 1/sqrt(C_mu), kappa² = (C_epsilon2-C_epsilon1) sigma_epsilon/Q.
For high-Re CLS, Q solves C(Q)Q²=1 and kappa solves

```
Q kappa^2 / sigma_epsilon + Y(R0) = C_epsilon2-C_epsilon1,
R0 = kappa Q^(3/2)/C_L,Yap.
```

In simple shear both closures give -R_xz = u_tau² and P_k = epsilon,
so the constant-k source is zero. The unforced momentum residual, using
the *retained* molecular viscosity, is

```
d_eta[(nu+nu_t) d_eta U] = -nu u_tau/(kappa r^2).
```

The dissipation residual (diffusive divergence plus implemented source) is

```
2 nu A/r^3 + [Y(R(eta))-Y(R0)] A^2/(K r^2),
R(eta) = R0 r/eta.
```

For standard k–epsilon Y is zero. For CLS the physical Yap distance remains
eta, hence the r/eta factor. These expressions follow directly from the
implemented high-Re equations: the transported CLS `eps_tilde` is used as
the high-Re dissipation variable and k sink. Selecting resolved `Molecular`
diffusion does **not** enable the low-Re gradient-k dissipation split or its
additional source terms.

Consequences:

- d=0 and zero resolved molecular diffusion recover the canonical bulk
  neutral equilibria. Finite polynomial approximation errors remain.
- A shifted standard profile is a bulk equilibrium with zero resolved nu,
  but the existing wall law/reference is not a displaced AWF implementation.
- Retaining nu leaves a momentum and epsilon residual even for d=0. The
  pure log profile is not an exact unforced finite-nu equilibrium. Its total
  matching traction also includes nu*d_eta U, whereas the reference AWF
  anchor is u_tau². A non-log total-stress solution or declared forcing would
  be a different reference.
- Displaced CLS generally has an additional Yap residual even with resolved
  nu=0. Replacing physical Yap distance by eta-d would change the model.
- Nonconstant k and finite-depth outer blends introduce further transport
  residuals. Exact preservation of the measured approach layer is not a
  consequence of selecting equilibrium AWF.

The native kernel probe checks 32 combinations of model, nu, d and height
against finite differences of the actual closure/diffusive fluxes. Its scaled
finite-difference tolerance is 2e-8; algebraic reconstruction checks use
5e-10. This derivative check is separate from solver preservation errors.

## Wall and experimental route

CV2 introduces **no displaced closure or wall-physics variant**. It retains
the current canonical models, physical wall geometry, physical strip inventory
0..hm, original momentum/production, Yap distance and integrated continuity.
The rough d=3.5 mm is an observational/data-fit parameter. `RANSHeightShift`
stays zero; its nonnegative additive log-law shift is not the paper's y-d.

The experimental route is a measured/inferred finite-depth inlet followed
by a computed empty-floor precursor. CV3 must calibrate its developed profile
at the would-be hill, with any canonical-model misfit retained. Molecular
resolved diffusion is the primary policy and HighRe is a paired sensitivity.

The smooth AWF uses ks=0, with its own momentum anchor
U_H=u_tau² J0(K,H,ks=0,nu). It must not inherit the rough-wall anchor. Both
closures retain the same EP derivation: mean production determines the capped
dissipation inventory, lambda=H*Pbar/A, cutoff=H*exp(1-lambda); the cap must
cover neither P nor H, wall Reynolds number must be admissible, and any
fallback is recorded. Original AWF is recovered exactly at zero weight.

## Versioned trial inlet and case generator

`cao_tamura.py generate` accepts a JSON configuration rather than changing
verification-pilot constants. It generates a provisional, unfitted case and
manifest, with a 0.8 m roof, 40 mm hill/100 mm half length or a paired floor,
physical normal H/P/N, explicit domain/mesh/span/resolution, canonical model,
viscosity policy, roughness and upstream transition. Every width uses the same
mesh. The narrowest planned ramp is resolved by 40 mm or smaller wall cells;
CV4/CV5 must establish sufficient resolution and a usable width plateau.

Only visible, positive-weight **Fig. 3a flat calibration data** initialize a
weighted cubic in log(eta/H). The outer mean blends with C2 weight to U_inf
between nominal 0.8 delta and 1.2 delta (0.2..0.3 m). k starts at Q*u_tau²;
the outer k=1.5*(I_outer U_inf)² assumes isotropy and I_outer=0.003 by default.
Epsilon is inferred from u_tau³/(kappa eta) and an explicit outer length scale.
Positive k and epsilon use geometric C2 blends. For original AWF an H..N
collar matches its epsilon trace; it is not a hill-data fit. All constants,
lengths and profile coefficients are recorded. The reserved Fig. 3d vertical
variance does not determine k. The unknown spanwise variance remains an
assumption to vary in CV3.

Re_H is based on U_H, not U_inf. The default analysis convention sets nominal inlet U_H=1
and nu=H_h/Re_H, then scales U_inf from the fitted flat U(H_h)/U_inf. This
is not a measured absolute wind speed or a verified developed incident Re_H.
CV3 must measure the developed U(Hh) and check that Reynolds number before
freezing calibration. An optional positive `assumed_nu`
provides a separately labelled physical dimensionalization. The rough trial
ks=6 mm is an unfitted input, **not** the cube height, z0 or displacement.

The top imposes U_inf, zero normal/span velocity, and nonzero outer turbulence.
The outlet uses zero prescribed residual traction through `RANSOutflow` and
natural scalar flux. These are declared starting choices, pending domain,
roof and facility sensitivity; the generator does not claim tunnel equivalence.
Its 3D periodic span is a proxy, not a representation of the end plates.

`prepare` constructs a native field from the actual named experimental inlet.
It uses the existing physical-normal/roof-compatible streamfunction, with a
Gauss-integrated inlet primitive in its N..2N blending band. Independent SciPy
quadrature verifies the input velocity/scalars, raw wall-normal condition,
projected positivity and span invariance. It rejects an inadequately resolved
positive inlet that projects to negative k/epsilon. No clipping is used.
The native preparer also checks finite positive geometry Jacobians.
Preparation tags the field as initialization only. The solver removes that
tag from its outputs; statistics also reject zero-time startup checkpoints.
A prepared field cannot enter statistics/calibration.

## Native physical observation

`RANSObservationAudit session.xml field.fld points.csv output.csv` reads the
native field with its metadata and exact mesh/expansions. Point input is
`id,x,span,z` in laboratory axes, including in 2D (span must be zero). It
supports tensor-product quad/full-hex meshes and their segment/quad boundaries.
It searches native elements, refines the inverse geometry map, and calls
`StdPhysEvaluate`; it never snaps to GLL nodes. Deterministic lowest global
geometry ID, then rank, selects a one-sided trace at shared faces. Outside
queries are reported and never extrapolated. Curved embedded edges use an
ambient-coordinate Newton inverse because Nektar's Geometry1D inverse only
implements regular segments.

The experimental floor also exposed nonfinite inverse seeds on nearly affine
curved quads. The production N sampler and observer now fall back to the
element centre for a nonfinite native seed and explicitly reject nonfinite
residuals. Previously, a NaN could pass a maximum-residual comparison as zero.
This fixes point location; it changes no closure or strip-continuity equation.

Resolved output uses laboratory velocity and Reynolds tensor components:
X streamwise, Y spanwise, Z vertical. In 2D native `v` is laboratory `w`.
Native saved model stresses are interpolated when present; otherwise the
instantaneous closure is evaluated and projected into the native expansion.
Molecular viscous stress is not included in Reynolds stress. Pressure is the
solver's modified pressure pi=p/rho+2k/3; AWF forcing reconstruction uses
mechanical pi-2k/3 and instantaneous convection.

Below H, the tool reads the same cached AWF K_P, forcing, roughness, continuity
derivatives and weight/cutoff as the production adapter. It reconstructs
the original auxiliary tangential/span momentum profile at physical eta.
The normal component uses the existing control-volume flux assumption,

```
J(eta) W(eta) = -eta d_s U_P - V(eta) d_span V_P,
J=1-c eta,  V=eta-c eta^2/2.
```

This is the implemented integrated-strip condition. It is not an assertion
of pointwise incompressibility of the separate locally parallel momentum
profile. The eliminated CV state U_P and auxiliary momentum velocity U_A(P)
are different quantities and have separate columns. The limiting normal
velocity difference from the resolved H trace is also reported. Tangential
matching is algebraic. The observer is read-only and does not open/truncate
solver wall-diagnostic files. The Python wrapper uses a separate session and
partition directory for each observation.

The auxiliary matching shear is checked independently with a fourth-order
backward derivative at H: D_A(H)*U_A'(H)=tau_A(0)+F*max(H-ks,0), for both
tangential components, both closures, smooth/rough walls and nonzero forcing.
The unforced streamwise anchor gives tau_A=u_tau², matching the resolved
HighRe neutral stress. Equality of the auxiliary and resolved *velocity
gradients* is not imposed: their eddy-viscosity laws differ. The finite-nu
resolved pure-log traction has the additional molecular term derived above.
For non-equilibrium curved stages, retain the measured weak tangential-flux
residual and refine it; the auxiliary ODE identity does not establish exact
pointwise shear matching of an arbitrary resolved solution.

Unavailable observations remain explicit:

- Below H, no full Reynolds tensor, k or epsilon profile is assigned. AWF
  supplies integrated budgets and selected CV/trace values, not a unique
  experimentally comparable full tensor. Those CSV entries stay blank.
- The compatible standard log wall has no below-H auxiliary profile here.
- Exact compact C1 feet have ambiguous one-sided continuity targets; their
  below-H probes are `unsupported_C1_foot`, not averaged into a made-up value.
- Auxiliary forcing is currently reconstructed for unforced Convective,
  non-dealiased cases. Body-forced/other-advection/dealiased auxiliary queries
  are explicitly unsupported; resolved native observations still work.
- Unknown experimental delta_75 prevents physical placement of Fig. 11.
  Dimensional near-floor amplitudes need an independent scale. Full scalar
  and ratio comparisons require the corresponding matched-flat/dividing-line
  sampling, not just the primary mean-profile probes.

Statistics use trapezoidal physical-time weights, with
R_total=mean(R_model)+covariance(U_resolved); no rms addition or diagonal
clipping. Sampling a transient does not establish stationarity. The `calibrate`
command evaluates a **flat-only candidate objective**, retaining all frozen
weights and unsupported coverage. It is a usable objective for CV3 parameter
search, not an assertion of completed fitting. `compare` evaluates supported
held-out profile families, optionally with same-ID matched-flat observations.
It rejects provisional manifests: CV3 must mark a completed, versioned frozen
calibration before reserved flat or hill scores can be inspected. Conditional
subset scores are never substituted for a complete-range score.
Statistics require consistent observer receipts/hashes and reject zero-time
startup fields as well as explicitly prepared/manufactured inputs. Solver
terminal L2 norms are not experimental comparison scores.

The speedup operator uses the matching surface/inlet/model/wall floor at x=0,
at the same **vertical** yprime. It checks the paired manifests, sampling window
and coordinates, and flags near-zero denominators. `probes --scope reference`
creates those same-ID floor probes; `local-reference` exports the additional
x-local floor data for drift assessment. `hill`/`reference` also include the
scalar crest point yprime/Hh=0.25. Longitudinal-rms and shear ratios use the
corresponding total covariances at matching heights, not rms addition.

`probes --scope dividing` supplies vertical near-wake traverses at x=0.11..0.40 m
in 10 mm steps. Starting tail heights are 2/4/6 mm for S and 7/9/11 mm for R;
`--tail-heights` declares alternatives. These are analysis choices, not claimed
experimental dividing-line heights. At each height the operator finds the
unique negative-to-positive mean-U recovery, without bridging missing probes.
It extrapolates x(y) from the lowest two complete tail heights to physical
ground and repeats with three for method sensitivity. Multiple recoveries,
missing support or fewer than three complete low tail heights are unavailable.
The output retains each finite-probe crossing, the ground-extrapolated x_R/Hh
and their distinction from a wall-shear zero or bubble length. An independent
synthetic line with x_R=0.237 m, a curved-tail perturbation and missing probes
verify the method without using the paper's x_R as a tuning target. Tail/grid,
height and sampling uncertainty must still be assessed in CV4/CV6. A separately
averaged native wall-traction trace is needed for the additional wall-shear-zero
diagnostic; it is not inferred from the dividing-line estimate.

## Restart and provenance

An experimental case has `RANSExperimentalSignature` and the expression
function `RANSExperimentalInlet`. The restart signature binds the declared
case fingerprint, actual named inlet expressions and referenced parameters,
and all `RANSExperiment*` observational parameters. Region 3 must prescribe
the same inlet expressions. This catches changed displacement/profile inputs,
in addition to existing geometry/roughness/physics/timestep history checks.
InitialConditions file paths and step count are excluded so full history can
continue across ranks. Legacy cases retain their original restart signature.

New field metadata `RANSObservationConfiguration` binds model/wall/geometry,
RANS parameters, viscosity policy and the experimental inlet. Observation
rejects a mismatched configuration or an old checkpoint without this binding.
Use the exact producing XML/manifest/mesh; this configuration fingerprint is
not a replacement for the recorded SHA256 of the producing session. Case,
data snapshot, binary and geometry-library hashes and actual process commands
are retained in workflow/validation receipts. Large native fields/logs stay
in the isolated Linux run root.

The native observer writes the verified configuration in a sidecar. Statistical
receipts require that configuration, experimental fingerprint, probe coordinates
and binaries to match. Producing XML filenames/step counts may differ between
compatible full-history segments; their individual hashes remain recorded.
Calibration/comparison checks that the resulting statistics belong to the
candidate manifest's physics fingerprint.
