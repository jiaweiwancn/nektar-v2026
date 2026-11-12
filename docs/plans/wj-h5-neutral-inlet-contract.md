# WJ-H5: neutral inlet and matching-height feasibility, zero timesteps

Baseline fetched unchanged: `e2fd85f0ef5d7587db7985ffc6670ae80df9992e` on
`codex/wj-earsm-awf`. Preserve all H2-H4 states, failures, ledgers and denials.
Old measured-rms/Blackadar/H_w=20mm validation, H4 capacity/refinement work,
long continuation, EP and physical hill runs stay paused. No governing step,
restart continuation, constitutive change or native capacity change is allowed.

New family: **WJNeutralBudgetReferenceV1**, an idealized full-height equilibrium
control, separate from WJFiniteDepthCompletionV2. Use exact Ne=81/20,
qe=81/(2sqrt(143)),Cd=572/6561 and
kappa=sqrt(1.30*(1.82-1.44)*sqrt(Cd)). Keep all original WJ constants,
full3D stress, DevStrainRawWorkV1, raw work and independent Cd scalar diffusion.
With H_h=40mm, U_infinity=1, use u_tau=.051, d0=.0875, delta=6.25,
ceiling20 and nu=.051*.005/3.97. The reported .2mm determines nu; fitted z0_e
does not change it. Molecular diffusion remains absent only in the outer model.

Original/worktree specifications must retain main SHA256
`1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124` and
appendix `0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c`.
Keep the quad/\quad typography erratum and correct (d+a*h)/a^2 coefficient.
The actual located paper is `/Users/jay/Documents/nektar-v2026/docs/literature/Experimental study on roughness effects on turbulent boundary layer flow over a two-dimensional steep hill.pdf`,
SHA256 `adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df`.
Pages6-7 and13-15 are read/rendered; do not publish the full source PDF.

Recompute physical heights from frozen Fig.3 CSV raw heights: panel a uses mm,
b uses z/delta, c/d use (z-d0)/delta. Keep reading resolutions as deterministic
reading bounds, not confidence intervals. Compute separate Fig.3 conservative
stress ceiling and the Fig.12 physical-z=7mm requirement (r=3.5mm). The latter
stays in every observation-region table even when it is inside the AWF.

Fit **one** outer intercept with fixed slope, u_tau and displacement, equal
weights in U/U_infinity, original physical-z centers <=25mm. Use
log(z0_e)=mean(log(z_i-d0)-kappa*U_i/u_tau) in consistent units. Frozen
sensitivities: omit the lowest of those three points; extend to centers<=50mm
(six points, excluding50.2788mm). Retain per-point residuals. Propagate the
declared height/value reading bounds separately and jointly, with membership
fixed by the original centers. Do not assign uncertainty to unquantified
d0/delta/u_tau anchors, tune windows or fit hill data.

Outer U=(u_tau/kappa)log(r/z0_e), k=K=qe*u_tau^2, epsilon=A/r,
A=u_tau^3/kappa. No measured-rms k or Blackadar/finite-depth tail enters it.
Record full covariance and the finite-depth mismatch at physical250/800mm.
At each H derive d_nu=a*H*z0_e/(H-z0_e), a=kappa*u_tau,
z_v=(nu-d_nu)/a and invert the unchanged negative-intercept roughness relation.
Keep printed z0=.2mm as a separate comparison branch. Reject absent inverses,
wrong branches and h_s>=H; keep inadmissible required-h witnesses in tables.

The planned candidates are H_w=5,5.5,6mm (dimensionless .125,.1375,.15), with
5.5mm central for organization only. Compute strict roughness/observation
intervals, 3% and2% outer molecular-ratio screens, H+, h_s/H and physical
plane/cube ratios. Neither reading bounds nor the 3% engineering screen proves
roughness-sublayer validity. The all-7mm-outside interval is expected empty;
no origin, h_s<H rule or observation inventory may change to avoid that conflict.

Inner U=u_tau^2*J0, D=nu+nu_t,A, P_ref=u_tau^4*nu_t,A/D^2;
Lambda=u_tau^4*(J0-nu*J2)/A, rc=H*exp(1-Lambda), Ec=A/max(r,rc).
Use nonconstant kw=K-integral_r^H Q/D, Q=integral_0^r(Ec-P_ref).
Independently integrate on physical breakpoints. Check endpoints, nonconstant
weak modes and full covariance. d_nu<2nu makes the appendix's convenient
proposition inapplicable; verify the unimodal-production extension conditions
P_ref(0)>A/rc, positive tail defect and zero total budget, plus actual profiles.
No source cap, clipping, flooring or evolved-state restoration is introduced.

Matching includes common U/K, **total** molecular+turbulent shear, zero k flux,
one dependent outer epsilon value A(K)/H and pi_w-pi_o=-K/Ne. Gradients need
not agree: U_w'/U_o'=1-z0_e/H. Use pi_o=0, pi_w=-K/Ne and compatible complete
regional outlet traction. ReferenceTractionV1 supplies top shear and epsilon
diffusive flux; top normal traction remains the impermeability reaction, not
an additional Function/pressure constraint. No extra gauge or H pressure row.

Use Functions on the flat CompositeTunnelV1 map, domain[0,8]x[0,20], Q3/Q1,
q12/normal16, IMEX1 single sweep, eta0, PimpleVerifyOff, positive wall nu,
uncapped Budget, endpoint K and StopOnBound. The stored positive TimeStep=.01
is a checkpoint identity only: **NumSteps=0 in every H5 session**. Fix one
horizontal cell. Filter all inner knots against H; include0,H and applicable
cutoff/roughness knots. Outer physical-r knots use H, existing measured mean
heights above H, 12.5-d0 and20-d0, then the exact composite inverse. No mesh
or timestep refinement campaign is allowed. Keep the3000-state guard.

At most three planned native CP0s, one per candidate, plus one extra ONLY after
a concrete generator fix. A model/representation failure is retained; no
unchanged retry or parameter/mesh tuning follows. Existing strict passive
readers may evaluate startup spatial residuals/weak moments, never solve or
invent a completed stage. Independent24/32 loads use original component
2e-9+2e-8*abs(value) and regional4e-9+2e-8*abs(sum) tolerances. Keep expression
approximation, analytical identities and native Q3/Q1 projection errors separate.

Linux volume UUID5AB8A669B8A642FD is verified at `/media/jay/新加卷1`.
New isolated root: `/media/jay/新加卷1/nektar-v2026-runs/wj-h5-equilibrium-inlet-20261112`.
Reuse governing producer c344dc2209f030e426e2ee5501d5e5b72a6d2c46 and native
observer e62b5f4d265f068dd5af3f090bfb0e385c3ae53e after source/runtime checks.
No native rebuild is expected. Preserve -ffp-contract=off -fno-fast-math.
One rank, all six numerical-library thread controls1; inspect current total
workers and remain below72. Keep artifacts outside Git except compact evidence.

New allocation:1200s total, including40s initial setup allowance and180s closure
reserve;8GiB process-tree RSS and8GiB new output. Reserve before launches;
unknown/failed work stays charged. Initial caps: focused algebra/calibration/
plots120s, each CP060s, each combined passive/native view150s, remaining source/
metadata80s. Reconcile measured costs within the fixed ceiling; no limit increase.
Reuse the existing bounded resource helpers; do not build a new provenance
framework or append to old ledgers. Zero governing attempts is an invariant.

Initial tolerances: exact constants/tensor identities relative~1e-12;
split reference integrals/weak identities scaled~2e-11 (absolute floors near
zero); expression approximation2e-11; native independent loads use the unchanged
mixed tolerances above. Rounded tabulated expectations use tolerances matching
their printed precision. Nonpolynomial represented residuals are measured, not
required to vanish or erased by forcing. Output explicit continuum/candidate/
measurement-scope/native-representation statuses and a bounded next decision.
