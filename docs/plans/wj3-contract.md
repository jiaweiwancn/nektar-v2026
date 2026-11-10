# WJ-3 native spatial and short-time comparison contract

This package continues reviewed `1d397855465be1daeefd160fc52c89c806ca0681`.
It does not replace the WJ-1 kernels or WJ-2 native single sweep. Mac owns
source/Git; isolated Linux builds and runs supply all numerical conclusions.
The two specification SHA-256 values remain
`1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124` and
`0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c`.
The literal `quad` in main eq:integrals-linear is the documented `\quad`
typographic erratum; specification bytes remain unchanged.

## Equations and selected route

The full three-dimensional WJ covariance uses the existing cubic physical N
root and DevStrainRawWorkV1 completion. Native 2D coordinates are physical x,z;
Ryy is retained. Praw=-R:Graw=Pdev-2k div(U)/3, and scalar diffusion uses
nu_d=Cd k²/epsilon independently of nu_eff. Defaults remain c1=1.8, c2=5/9,
Ce1=1.44, Ce2=1.82, sigma_k=1, sigma_epsilon=1.3, Cd=C*=572/6561,
qe=81/(2 sqrt(143)), kappa=.381917863546733. Outer molecular diffusion is zero;
wall molecular nu=1.5e-5 is positive throughout.

The inner stress remains full-strain Linear. EP has uniform chi=1, endpoint
K=k(H), zm=H/2 and complete uncapped source P-pA0. Budget has zK=zm=H,
source P-Ec and the transported nonconstant inner k reference. Reduced pressure
is pi=p/rho+2k/3; its single global gauge preserves the regional mechanical
jump -K/Ne. Shared velocity/k, distinct gradients and pressures, opposing
normals and one dependent outer epsilon value trace remain unchanged.

IMEX1/PimpleSingleSweep performs one vector-U, one pressure/gauge, one k,
one outer-epsilon solve and one logical MassProjection per completed attempt.
There is no equilibrium forcing, anchoring, clipping, reset, condensation or
nonlinear correction. eta=0, gravity/body/scalar forces=0, EvaluationReuse=Off,
PimpleVerify=Off. WJ now rejects Verify=On before integration because that
legacy option performs additional governing verification solves. All prior
unsupported model/search/transition/cap/curved/3D/MPI routes stay unsupported.

## Immutable measurements

`wj_measure` extends the strict checkpoint reader. Its entry alpha is exactly
zero and time equals the checkpoint time. Every GLM slot/time label, primary
coefficient and external boundary value is verified. It does not form or advance
an integration callback. The solver report is SHA-bound and must match the
checkpoint state, physical signature, time and step. For an evolved state its
last_stage_history/time/alpha are additionally matched to the accepted WJ audit
packet's hashes and context. A startup has no completed-stage residual.

Measurements retain separately:

1. Uncancelled physical spatial residual at alpha_mass=0, eta=0.
2. Completed-stage residual with the actual last implicit callback history.
3. Drift relative to this mesh's unchanged represented initial field.
4. Total error relative to the correct piecewise continuum profile.
5. Same-time differences on a common physical observation rule.

Residual norms use the existing consistent native mass dual and fixed physical
scales. Their blocks include u,v,k, the entire outer epsilon region, inner and
outer continuity, gauge and boundary/trace constraints. Source, diffusion,
interface and regional load vectors are retained; their separately computed
norms cannot be added as if they were signed vectors. Synthetic WJNativeChecks
history=state, alpha=1e6 remains a separate frozen-assembly test only.

Physical L2 and sampled Linf errors use one common rule split at the finest
supported element boundaries and zc,zv,h,zm,H where applicable. Both H sides
are evaluated. A single volume-weighted pressure gauge shift applies to both
regions and to mechanical pressure. Unequal coefficient vectors are never
subtracted. The initial projection's prescribed interface derivatives partly
impose startup traction/flux agreement; these initial identities are labelled
separately from evolved measurements. A discrete first-step pressure adjustment
is allowed; its mesh dependence and common-time convergence are measured.

## Audit witnesses and arithmetic observations

WJFiniteAuditV2 reports separate scalar-domain, profile-gate/budget, covariance,
root-branch/residual, scalar-viscosity, raw-work and two principal-form conditions,
with individual margins and failure categories. Each limiting principal witness
includes NativeElementSideV2 geometry/element/side, physical point, k, epsilon,
raw G, sampled K, N, covariance spectrum, nu_eff, nu_d, direction xi, polarization
a, a·xi, full dR and H, contraction and normalization. State/history hashes,
phase, time and alpha bind the witness to its observed stage.

The in-plane search uses 32 periodic brackets and 48 golden iterations per local
minimum. Transverse polarization uses the symmetric 2x2 eigenvalue for x-z
propagation. Independent complex differentiation and bounded angular replay
check limiting witnesses. These finite spatial/plane checks and the exact
reference nu_eff/250 certificate are not general-state/global stability guards.
Both adverse WJ-1 states remain in the audit and regression suite.

Optional PimpleArithmetic copies already assembled matrices, right-hand sides,
solutions, residuals and mass/correction maps at a requested physical time. It
performs no factor application or governing solve. An On/Off trajectory check
requires bitwise equality. Offline SVD inspects conditioning; long-double
matvecs estimate linear backward errors. The declared comparison-floor estimate
combines these conditioned errors, n*epsilon arithmetic formation reserve,
physical-mass conversion, pressure-to-velocity correction propagation and
measured alternative field-evaluation sensitivity. It accumulates 2*N terminal
step estimates and applies safety factor 100 before interpreting rates.
This is a finite short-trajectory numerical estimate, not a rigorous nonlinear
perturbation bound. Spatial residuals and nonlinear stage defects are excluded
from the arithmetic floor. Unresolved differences receive no observed order.

## Bounded execution design

Fresh startup data are projected independently on p=3 meshes formed by 0/1/2
bisections of the unchanged graded normal intervals (expected 639/1257/2493
unknowns). The dense 3000-unknown limit remains. Source quadrature 16 versus32
and independent observation orders6 versus10 are checked separately; basis
quadrature stays8x8. Independent source/load integration splits physical and
element breakpoints. Rough-reference hypotheses zv<=0, d>=2nu, Lambda>=1 are
checked separately from represented-field acceptance.

Scouts use dt=1e-6,4e-6,1.6e-5,6.4e-5,1.28e-4, up to four fresh steps each,
stopping at the first hard failure or material projection. Each dt owns a fresh
native CP0 and native signature; primary/gauge and corresponding physical
histories must agree. No signature or history is edited. After cost measurement,
the driver writes a frozen D,T=32D design before any 32/64/128-step main run.
It records a practical 3600-second per-policy remaining trajectory budget.
The same-T cross-mesh run uses2493 if a conservative cubic-cost estimate fits;
otherwise639, leaving finer-mesh temporal qualification open. The total cap is
500 scheduled attempts per policy, including failed attempts. Work is serial,
one build worker, one test, numerical-library threads=1, below the server-wide
72-worker ceiling. No historical long campaign is included.

Two common-time differences A and B produce p=log2(A/B) only above their
declared arithmetic floor. Total continuum error, initial drift, spatial
residual and completed-stage defect are distinct outputs. A hard obstruction
retains its accepted predecessor, failed candidate and complete witness; the
frozen target is not quietly shortened or repeatedly halved. A 16+16 restart
must equal the 32-step primary/gauge/history/time, accepted audit and samples
exactly. Later frozen readback and end-of-step rejection check unchanged saved
bytes, rollback and honest solve counts.

These gates establish measured spatial consistency and resolved short-time
behavior where they pass. They do not establish exact finite-space equilibrium,
long-time preservation, global positivity/realizability, coupled nonlinear
stability, native 3D/MPI or terrain qualification.
