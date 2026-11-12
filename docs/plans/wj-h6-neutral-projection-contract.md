# WJ-H6: neutral-reference projection diagnosis, zero evolution

Fetched starting HEAD and origin both equal the reviewed
`34937bb2ed8730bd836ff63f13e278d3e3facada`. The WJ worktree is clean;
unrelated worktrees and the main checkout's untracked literature are preserved.
No published history or previous run root is rewritten.

H6 diagnoses H5 initialization and executes exactly three fresh CP0 launches
on one prescribed mesh. All sessions have NumSteps=0,Time=0; TimeStep=.01 is
stored metadata only. No restart, completed-stage packet, governing timestep,
extra CP0 retry, p-order/capacity/backend change or alternate mesh is allowed.
The old measured-rms/Blackadar/H=20mm work, H4 investigations, EP, long
continuations and physical hill runs remain paused.

Keep WJNeutralBudgetReferenceV1 with Ne=81/20,qe=81/(2sqrt(143)),Cd=572/6561,
kappa=sqrt(1.30*(1.82-1.44)*sqrt(Cd)), the original WJ coefficients, complete
3D covariance, DevStrainRawWorkV1, raw stress work and independent Cd diffusion.
Keep u_tau=.051,nu=.051*.005/3.97,d0=.0875,delta=6.25,top=20,H_h=40mm,
K=qe*u_tau^2. Viscosity still derives from reported .2mm/3.97, not fitted z0.
Molecular viscosity is omitted only in the outer operator and retained in AWF.

Reuse the exact H5 fitted-intercept binary64 value (display
.227050812848672mm, dimensionless .0056762703212168), degree16 XML Functions,
finite-nu nonconstant inner k, current-K counterfactual cutoff, dependent A(K)/H
epsilon trace and pressure offset. Keep the compatible complete outlet
tractions, ReferenceTractionV1 top shear/epsilon flux, physical coordinates and
one pressure level. No extra gauge, inner epsilon PDE or H flux is added.
H=5,5.5,6mm and h_s=4.381467178582103,4.359758967365407,4.341801417014184mm
remain unchanged;5.5mm is an organizational center only.

Retain every observation-region entry. Fig.3 requires H<6.298704214mm and
h_s<H requires H>4.413515542mm. The all-7mm-outside requirement H<3.5mm has
no solution. Both reading intervals stay distinct: bounded refitted-estimator
range[.2171963971,.2370377766]mm and intersection of individual point-compatible
intervals[.2261849783,.2309394239]mm. Neither is a confidence interval.
The finite-depth mismatch and roughness-sublayer assumption remain unresolved.

The original/worktree specification hashes are unchanged:
main1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124,
appendix0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c.
Retain the literal quad/\quad erratum and correct(d+a*h)/a^2 coefficient.

Independent replay stages are: exact continuum/XML evaluation; physical
unconstrained L2 projection with the native mapped order16 mass/load rule;
OpenConstraints exterior coefficient lifting; actual represented-K dependent
epsilon trace lifting; retained final native coefficients. The middle stages
are **reconstructed**, not native-exported. Boundary transforms use native
Modified_A basis,8-point GLL endpoint-constrained trace equations and exported
ids/signs. Pressure remains regional Q1 and unconstrained by Dirichlet data.
The initial default order8 canonical projection followed by physical corrections
is checked through its final physical order16 mass/load equations, without
claiming to export native iterative intermediate values.

Attribution is permitted only after reconstructed final coefficients/owned
fields/gradients agree with native at scaled2e-11. Field scales are U=1,
k=K,epsilon=A/H,pressure=u_tau^2; coefficient comparisons have1e-15 absolute
floor. Physical velocity gradients use U_infinity/H_h=1 with2e-11 absolute
floor. Local principal comparison uses2e-11 absolute tolerance, root residual
1e-12. Record linear backward residuals and XML-versus-continuum errors
separately; expression accuracy remains2e-11. Never divide by a zero field.
Passive G-only and k/epsilon-only replacements diagnose sensitivity and are
never saved as initialized states. Rejected H5 candidates retain their native
report/event role and are not admitted to the checkpoint reader.

The sole new mesh bisects every original chart interval once, including both
regions. Preserve all old endpoints, H,cutoff/roughness and observation faces;
verify physical midpoint equivalence under the affine regional maps. Nx=1 on
[0,8],Q3/Q1 and native points remain unchanged. Expected Ni/No/N are8/17/1174
for H5 and16/34/2328 for H6, using N=38Ni+50No+20. Verify these before launch
and from actual native output. The3000-state guard remains in force.
All H5 native settings remain fixed: q12/normal16,Functions/Physical,
CompositeTunnelV1,OpenTunnelV1,ReferenceTractionV1,IMEX1 single sweep,eta0,
PimpleVerifyOff,uncapped Budget,InterfaceTrace,StopOnBound and the same reuse
setting. H5 generator defaults and historical Function bytes must reproduce.

Freeze the common physical comparison inventory by rule, per height: union
of baseline/refined element faces, outer normal Gauss12 and native LGL8
locations at x=0,4,8; H,all retained Fig.3/Fig.12 points; all exported H5/H6
principal witnesses;1e-6 physical inward/normal offsets at those witnesses
and H; and exported boundary/interface audit rows. Reconstruct every
legitimate owner at a shared physical point. At H evaluate WJ/epsilon only
on the outer owner and AWF stress/local k on the inner owner. This finite
comparison inventory is separate from each run's full native hard inventory.
The inventory rule is frozen before launch; final witness coordinates merely
instantiate it, without selecting another mesh or discarding an earlier point.

Gates remain separate. A: unchanged reference/coordinates/Functions/finite-nu
conventions. B: all original native physical predicates. C: every declared
outer common-point owner has in-plane m/nu_eff>=.5*m_exact, with m_exact
computed by the independent action (~.004137527847); positive transverse and
covariance margins separately. Define E_G=max|Gxz/Gxz_exact-1| and
E_epsilon=max|epsilon/epsilon_exact-1| over each shared outer inventory.
Both refined values must decrease versus H5. One bisection is sensitivity/
recovery evidence, not an order or global ellipticity proof.
D: accepted refined checkpoints pass strict read-only observations, exact
shared/dependent value traces, and q12 versus independent24/32 weak checks
using actual projected native test traces. Keep component2e-9+2e-8|value|
and regional4e-9+2e-8|sum| limits. Report physical spatial mass-dual residuals,
traction/k-flux defects and constraints without forcing zero. H5 rejected
startup norm packets are not comparable accepted-reader mass-dual norms;
mark unavailable comparisons explicitly. Executed unresolved weak checks are
UNRESOLVED, never NOT_RUN. No completed-stage context exists in H6.

Linux authoritative root, after verifying UUID5AB8A669B8A642FD:
`/media/jay/新加卷1/nektar-v2026-runs/wj-h6-neutral-projection-20261112`.
Reuse governing producer c344dc2209f030e426e2ee5501d5e5b72a6d2c46 and
observer e62b5f4d265f068dd5af3f090bfb0e385c3ae53e after source,ELF and loaded
library checks; no C++ change/build is expected. Preserve RANS_BUILD_SOLVER=ON,
-ffp-contract=off,-fno-fast-math. Use one rank and OMP/OpenBLAS/MKL/BLIS/
NumExpr/VECLIB threads1; inspect current capacity and stay below72 workers.

Fresh persistent allocation:1200s including40s setup and180s reserved closure;
8GiB process-tree RSS and8GiB new output. Reserve work before launch, retain
failed/unknown reservations, and never append to H5 accounting. Caps are120s
for each focused/reconstruction operation,90s per CP0,150s per combined passive
view and80s publication. Exactly three unique-height CP0 reservations, no
retry exception. Continue remaining planned heights after an individual failure
only if resource gates allow. Targeted tests include the nine H5 methods and
H6 bisection/counts,polynomial mass replay,lifting/actual-K trace,lineage,
projected-exterior-test and honest status classification. No broad replay.

Deliver compact stage/quality/weak/interface/accounting tables and scientific
plots. Preserve heavy artifacts and all intermediate reconstructions in this
new root. Close PASS/OBSTRUCTED/UNRESOLVED/NOT_RUN separately. A successful
startup study can propose a later bounded neutral-control evolution package;
H6 itself executes no such extension.
