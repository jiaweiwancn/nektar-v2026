# WJ-H1: budget-only composite rough-hill foundation

Frozen before native scientific execution, 11 November 2026. Reviewed local
baseline is `b47fd90ce5a67b5afd9af8378ad1f139976b98b3`; the same GitHub head
was independently read through the Linux connection. Mac HTTPS fetch attempts
timed out; fetch/push receipts will retain the actual transport used. No old
WJ-4/5/6/7 evidence, ledger, checkpoint or runtime is changed. WJ-7's joint
RESOURCE_LIMIT (Budget 10866.759210 > 10800 s), zero new steps and T3 NOT_RUN
remain closed. Hill states have a new physical signature and fresh histories.

## Case and interpretation

The experimental target is Cao & Tamura, JWEIA 94 (2006) 1–19,
doi:10.1016/j.jweia.2005.10.001, rough hill in rough inflow and a separate
rough flat reference. A local original PDF and extracted Markdown exist in
the original Mac literature directory. A machine-readable case records page/
figure provenance, quoted rounded anchors, conversions and assumptions. No
figure profile is invented. The finite-depth Fig.3 inlet, missing spanwise
variance/epsilon completion and physical inlet/outlet qualification are later
work. The exact neutral WJ variances are not retuned to the measured rms.

Scale lengths by H_h=.040 m and velocities by U_infinity. Hill height is1,
half-width2.5; physical cube height .125, spacing .875, density .041,
aerodynamic z0=.005, assumed displacement d0=.0875, delta~6.25, ceiling20,
u_tau=.051 and positive wall nu=.051*.005/3.97. Cube height, AWF sand grain,
aerodynamic intercept and displacement are distinct inputs.

The *provisional verification* candidate is constant normal H_w=.375 (15mm).
Invert the existing AWF correlation on a flat wall using the WJ logarithmic
slope and rounded z0/Re_s, not hill speed-up or reattachment. Its expected
h_s~.090675 (3.627mm) is conditional. Report inversion, at least two candidate
matching heights, cutoff/gate, positivity and molecular/turbulent ratios.
No experimental inlet calibration is claimed. A controlled small-velocity,
positive-k/epsilon, reference-Dirichlet verification problem is permitted;
its fields/forcing are labelled separately from an experimental inflow.

Future physical boundaries: finite-depth inlet; downstream pressure/traction
level with scalar outflow and specified backflow; ceiling approximation at
the original physical Z=20 (slip is a candidate needing sensitivity). A
pressure-level outlet replaces, rather than supplements, a global gauge.
The present closed Dirichlet prototype retains one global pressure gauge.

## Geometry and operators

Keep explicit faces at s=+-2.5, including both one-sided curvatures. Let
f(s)=cos^2(pi*s/5) inside those feet and zero outside, g=sqrt(1+f'^2),
n=(-f',0,1)/g and kappa_phys=-f''/g^3. The effective ground is
X_eff=(s,0,f)+d0*n; it is not a vertical displacement on slopes.
Inner X(s,r)=(s,0,f)+(d0+r)*n, 0<=r<=H_w. Relative effective-wall
curvature is kappa_eff=kappa_phys/(1+d0*kappa_phys); base area is
g*(1+d0*kappa_phys), and J=1+r*kappa_eff. Native physical weights include
this geometry exactly once. Source ray weights use base area times J dr.

Outer canonical rho in [H_w,20] uses eta=(rho-H_w)/(20-H_w) and
Phi=(1-eta)*X_H+eta*(s,0,20). Eta is not normal distance. At H the physical
interface and tangential derivative coincide, but normal computational
derivatives differ; physical gradient evaluation owns its native element.

The normal tube is certified only through d0+H_w. With
B=pi^2/(2*2.5^2), require (d0+H_w)*B < .9. This retains the legacy .1
minimum stretch and global C1,1 injectivity argument. For a=1+(d0+H_w)*
kappa_phys, the outer determinant times (20-H_w) is
(1-eta)*a*[20-f-(d0+H_w)*g]+eta*[20-f-(d0+H_w)/g].
Certify both brackets positive and a>0 using whole-domain analytic bounds.
Horizontal derivatives stay positive; at fixed physical x the admissible
eta interval is connected and physical z increases, giving injectivity.
Validate the whole tunnel, simple boundary, interface/cofactors and both feet.
No full-height normal-ray extension or relaxed reach certificate is allowed.

For completeness the two region interiors do not overlap. The offset interface
is a graph with slope f' and global Lipschitz bound M. Each outer segment obeys
z(eta)-z_H(x(eta)) >= eta*[Z-z_H(s)-M*|s-x_H(s)|] > 0, because
(1+M*|f'|)/sqrt(1+f'^2)<=sqrt(1+M^2) and the certified ceiling bracket is
positive. The normal tube lies below that interface. This also certifies the
inclined-plane union when its lateral chart edges bend at H.

Keep full 3D WJ tensor/cubic surrogate, Original, DevStrainRawWorkV1,
physical epsilon, Praw=-R:Graw and independent Cd scalar diffusion. Constants
remain c1=1.8,c2=5/9,Ce1=1.44,Ce2=1.82,sigma_k=1,sigma_epsilon=1.3,
Cd=572/6561,qe=81/(2sqrt(143)),kappa=.381917863546733. Omit molecular
diffusion only in the outer high-Re operator. Inner full-strain Linear stress
uses local transported k, positive nu, endpoint K coefficients and source
P_A-E_c. The existing metric cutoff solves integral(J Ec)=integral(J p_ref)
at current K. No cap, chi, redistribution, inner epsilon PDE, extra epsilon
flux, clipping or reference reset. One dependent epsilon trace A(K)/H_w.

Q3 U/k/epsilon and regional continuous Q1 pressure retain shared U/k,
separate pressures/gradients and one gauge. Pi=p/rho+2k/3. Full physical
stress/traction and tangential projection of sampled outer mechanical
pressure/convection drive the inner roughness force. Its signed work is
reported. Normal pressure jump is a consequence of common weak traction,
not an extra pointwise pressure condition. Preserve the PIMPLE pressure
operator D*(alpha*M)^-1*G-C, consistent mass and actual roughness G.

## Bounded implementation and Linux gates

New heavy root: `/media/jay/新加卷/nektar-v2026-runs/wj-h1-20261111`.
Use isolated immutable source epochs/build, installed dependencies read-only,
RANS_BUILD_SOLVER=ON, Release -ffp-contract=off -fno-fast-math. One compile
worker, one serial native process/test, six numerical-library thread controls1;
inspect current workers/memory/storage before launches and stay <=72 total,
including application/MPI helpers and other campaigns. Build ceiling7200s;
numerical study ceiling10800s; <=40 scheduled prototype/control attempts,
<=3 accepted steps per terrain trajectory, native invocation<=900s,
RSS<=8GiB, invocation output<=4GiB, H1 output<=16GiB. Retain failed work and
actual time/RSS/output costs. No large SVD or all-history campaign.

The narrow H1 CTest runner reuses the existing process-tree/worker inspection
helpers but opens only `wj-h1-*/checks-accounting.json`. It reserves every
selection before launch, keeps failed/unknown reservations charged, and bounds
whole-process RSS, current worker count and output. A600s closure reserve and
120s setup allowance supplement actual check wall time. Native timing uses
the existing passive timing stream, with per-invocation source/runtime/argv
receipts. No old package ledger is opened or reset. Build time is separate.

Initial coarse chart mesh retains feet and crest, s=[-5,-2.5,0,2.5,5],
rho=[0,H_w/2,H_w,2,6.25,20], p=3; confirm dimensions before evolution and
keep the legacy3000 dense-state guard. Start with dt=1e-6 in H_h/U_infinity
units. Controlled flat, inclined and curved cases use the same declared
state scales; these are implementation checks, not developed hill predictions.
If a demonstrated representation/implementation failure requires a changed
prototype, retain it and record the named change before another launch.
Never retry an unchanged hard-rejected state or change the model to pass.

Required fresh checks: full-height geometry/zero hill/incline/physical-gradient
FD/volume/surface/cofactor identities, invalid geometry and both feet; existing
WJ/Budget kernel oracles and unchanged small flat native regression; signed
traction/reversal; curved metric budget and pressure action versus residual/
frozen columns; independent source weights/interface modes; 1–3 actual coupled
native steps, all native hard predicates and inactive MassProjection. Preserve
partial solve counts/candidates on failure. Report true pointwise traction/k
flux defects separately from conservative shared-test assembly. No zero-step
readback or kernel-only result counts as the terrain demonstration.

The inclined-plane control reverses the initial vortex sign. Before any native
governing launch, the first compiled geometry driver rejected its incorrectly
placed legacy AWF reporting point P=H. The correction uses P=H/2 only in that
scalar-kernel inversion; the physical WJ endpoint K and epsilon trace stay at
H. The first coupled-probe compile also lacked explicit diagnostic/hash
headers. Both failed source/log epochs are retained; neither is a model-domain
failure or an executed governing step. The verified installed CMake package is
`/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake`.

Record physical tensor axes and both interface owners, K/k/epsilon, cutoff and
metric margins, covariance/root/principal diagnostics, actual/reference
production, wall traction, global budgets and correction accounting. Point
observations use original physical coordinates; experimental y'=z-f(x) is
vertical clearance, not r. Crest speed-up .875 at y'=10mm, reattachment6.5H_h,
and wake U at physical z=7mm are future observations, never fit targets here.

Terminal components are PASS, OBSTRUCTED, UNRESOLVED, RESOURCE_LIMIT or NOT_RUN.
H1 success means a verified mapped coupled prototype only. Estimate future
hill/wake resolution, development clocks and dense bottlenecks without
claiming scalable capacity, calibration, experimental agreement, stationarity,
global positivity/realizability/stability, native3D/MPI or terrain convergence.

The controlled fresh initial velocity is the curl of
psi=.1*r^2*(H_w-r)^2/H_w^3 in the inner region and zero outside. In physical
components U=.2*t*(1-t)*(1-2t)*tangent, t=r/H_w. Initialize its shared H trace
to zero once; retain projected exterior differential traces thereafter. Inner
k=K0*[1+.1*(1-3t^2+2t^3)], outer k=K0 and outer epsilon=A(K0)/H_w.
This is positive verification data, not the budget stationary reference or a
finite-depth experimental inlet. It is never reset during evolution. Run three
hill and zero-hill steps and one inclined-plane step using the production
coupled patch and native IMEX1 history implementation. Keep full solver build
and an unchanged small flat native Budget regression as a separate control.
Independent evolved loads use three declared shared modes (1,s/5,an H hat),
with retained exterior reactions; full frozen blocks and affine pressure
columns receive separate native checks. No all-mode nonlinear derivative claim
is made for those three weak-load modes.

Specifications remain byte-identical to the originals:
main1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124;
appendix0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c.
Retain literal quad/\\quad erratum and correct (d+a*h)/a^2 coefficient.

The first flat Budget regression completed its native transaction but failed
the old curved-selector message assertion. The new named terrain route still
rejects that cylinder input; update only the expected exact reason. Retain the
failed CTest receipt and its five charged attempts. Re-run with the unchanged
source-v3 native binaries and a separate CTest metadata epoch pointing to the
corrected Python fixture; do not overwrite the original source/build metadata.
The additional five attempts keep the planned total at32, below40. All current
native sources remain at b50aaa8; later publisher/fixture commits do not require
recompiling or reopening the already-passed terrain trajectories.

Sealing found a separate loader defect in the inherited optional prebuilt
SpatialDomains path: CTest always requested a new-build geometry-lib directory,
which did not exist in that mode. The loader used the installed library rather
than the intended isolated source-built one. Preserve native-v3, summary-v1,
its default-shell ldd receipts and all results as the installed-library lineage.
Fix the execution path to the actual prebuilt library directory and register a
zero-solve runtime-closure check. Native C++/headers are unchanged; reuse their
verified ELFs, configure fresh test metadata and repeat only the three mapped
controls (3+3+1 steps). Their seven extra units make39/40. No repeated flat or
legacy transaction is authorized/needed for this CMake-only correction; retain
their completed checks explicitly under the original installed dependency
closure. Recheck relevant zero-solve oracles and a flat accepted-checkpoint
readback with the corrected loader. Record both dependency identities. New
mapped output/observer epochs must not replace their predecessors. The H1
numerical ceiling remains10800s and the default evolution size limit3000.
