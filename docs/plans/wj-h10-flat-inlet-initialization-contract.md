# WJ-H10: compatible flat-neutral startup and three first steps

Fetched baseline is d5c708457c5248b2bb5e719aebd72628e007a8b7; there are no newer
WJ commits or dirty WJ files. Main-checkout untracked files and the Linux shared
checkout's unrelated docs edit are preserved. No existing H10 roots/jobs were
found. Use new ext4 `/home/jay/build/wj-h10-flat-inlet-20261113`; verify bounded
write/fsync/readback, free space and live allocations before launches.

Implement only opt-in `RANSWallInitialVelocity=FlatInletTraceV1`; absent/Off keeps
legacy signatures and behavior. Require serial2D WJ Budget/OpenTunnel, physical
Functions, horizontal zero-amplitude Plane CompositeTunnel, constant H and
conforming native Q3 Modified_A tensor partitions. Reject inclined/curved maps,
nonzero essential normal velocity, inconsistent native partition/IDs and any
initial perturbation. Reuse the already lifted native inlet coefficients: local
radial coefficients include inlet signs; copy them to horizontal endpoint modes
0/1 with destination signs, zero higher horizontal modes and free normal velocity.
No analytic reevaluation, RadialProject, mass/KKT/pressure solve or new flow
assembly. Essential coefficients, including signed zeros, remain bitwise fixed.

Here physical x=chart x; z=d0+rho inside and
z=d0+H+(rho-H)*(Z-d0-H)/(Z-H) outside. The injected U_h(rho),W_h=0 therefore has
physical div U=0 and all regional Q1 moments zero in exact arithmetic. Shared
normal faces/H retain velocity continuity. Natural outlet velocity is free.
Apply after Project/Lift/SetReference and m_initial freeze, before integrator
history initialization. Fixed reference/defect, k/epsilon and pressures never
change. Keep old residual fields attached to m_initial; add separate actual
pre/post state diagnostics and invocation count. Restart/strict readers apply
zero transforms; mode/version enters patch/shared checkpoint identity only when
active. Default/explicit Off signatures are identical.

Freeze the H6-refined6mm family: H=.15,physical zH9.5mm,hs4.341801417014184mm,
d0=.0875,u_tau=.051,nu=.051*.005/3.97,16inner+34outer cells,Nx1[0,8],Q3/Q1,
2328 states,q12/normal16 and3000-state guard. Clone H8's hash-verified6mm XML:
physical Functions/geometry/constants/reference bytes remain exact. Keep original
full3D WJ/DevStrainRawWorkV1/raw P and Cd diffusion, wall nu/high-Re outer,
nonconstant inner k,current-K counterfactual cutoff/dependent epsilon, original
pressure seed, complete outlet/reference-top tractions/scalar fluxes, IMEX1,
PimpleSingleSweep,k_then_epsilon,eta0,ExactStepV1,PimpleVerifyOff,StopOnBound.
No equations, pressure space, forcing, clipping, retuning, capacity or mesh change.
Specs retain main1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124,
appendix0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c,
quad/\\quad erratum and correct(d+a*h)/a^2 coefficient.

Full-size matrix, in order .01,.005,.0025: exactly one fresh NumSteps0 CP0,
one strict CP0 view, only after CP0 qualification one restarted governing step,
one strict accepted CP1 view. Max3 CP0 launches,3 attempts,6 views; failed/unknown
reservations consumed, no retry/substitute. First .01 CP0 is retained cost pilot.
Require all CP0 physical/fixed-reference arrays bitwise equal with valid own dt
signatures/history labels. Each restart reports zero Project/Lift/initializer calls.
Enable existing first-time PIMPLE arithmetic and k snapshots on each governing
step; do not alter PIMPLE or add solves. Small native synthetic fixture is separate
from full-size cases and cannot become a CFD trajectory.

CP0 gates: all34inner+70outer signed Q1 moments<=2e-15; all declared owned physical
points |div U|<=1e-12*max(1,||G||F); exact essential/nonvelocity/reference/defect
preservation; unchanged current K/cutoff/dependent trace; original native scalar,
root,covariance,diffusion,ground/outlet/raw-work admission; zero bound activation,
fallback or material correction. Expanded principal>=.002068763923454 (half
exact .004137527846908), positive transverse/covariance separately. No requirement
to beat old CP0 .002544936550653 or attain exact continuum fields. Inventory is
H8's complete6mm frozen union, all face/H owners and established local witnesses;
add new native witnesses with the existing1e-6 physical neighborhoods. This finite
scope is not global ellipticity. Hard/quality exit or unresolved sign within
2e-11 halts further full-size launches, retaining original candidate/history.

For each accepted CP1 recheck all gates except initial exact-divergence (record
its evolved divergence), plus solve counts1U/pressure/k/epsilon/one bound CALL,
zero activations and raw/corrected k equality. Bind P,C,F,K,A/true callback,
alpha/time; expect P=C static lifting, field masks and F/A flow equality. Use
saved momentum/pressure matrices, vectors, row maps and BGp only for matvecs:
verify du=predictor-BGp*dpi, pressure_rhs=D*predictor+rp(C),
rp(F)-rp(C)=D*du and rp(F)~=0. Complete Gp includes roughness/exterior pressure.
No extra governing factor/solve. Diagnostic pressure-mass algebra is counted.
Compare signed dpi,j=dt*dpi,predictor,correction,du,predictor/dt with H9 without
regauging, full/regional physical norms/means and common-owner profiles. First
steps occur at different times: no order/common-time convergence claim or forced
pressure target. Old unsaved predictor/correction vectors remain unavailable.

CP0 and3 endpoints get existing projected1,x/5,interface_hat signed weak loads,
state-dependent reactions/current-K recomputed. Deduplicate identical CP0. Use
q24/q32 and unchanged component2e-9+2e-8|value|,regional4e-9+2e-8|sum| limits;
keep additive terms, aggregates, diagnostics, spatial and actual-stage defects
separate. No all-intermediate weak campaign or quadrature search.

New persistent allocation3600s: setup/build/tests/runtime<=2400 including fully
charged60s setup; science<=900; fully reserved/charged300s closure/transport.
Reserve before operations, keep failures/unknowns/monotone child halts. Initial
serial build estimate1489s; cap/remaining admission uses25% headroom plus full
closure. At first .01 CP0 update from actual build/pilot cost. No cap increase,
reset, borrowing or uncharged retry. Limits4GiB tree RSS,8GiB isolated build/source
footprint,2GiB science output checked separately; single compiler --parallel1,
serial numerical worker/all six numerical thread controls1,live total<=72.
Newly built governor/reader have new loader/library/ELF identities; qualify
prebuilt SpatialDomains/dependency ABI. Never attach old H8 ELF identity to new
binaries or mix private glibc closures. Kernel reuse only with its verified wrapper.

Run51 existing H5-H9 methods plus focused coefficient/sign/map/partition,
preservation/idempotence, default/signature/restart/history, arithmetic/budget
fixtures; report actual execution. Preserve old H7/H8/H9 ledgers, denials and
UNRESOLVED scientific conclusions. Commit/push normally and verify remote HEAD.
Deliver compact versioned evidence; heavy matrices/fields stay in this root.
If all gates pass, recommend only separately bounded continuation to commonT=.04;
none is executed here. Finite-depth inlet/measurement scope, physical H selection,
7mm inside AWF/empty all-outside interval, hill, EP, old measured-rms/20mm,
capacity/backend and longer work remain pending.
