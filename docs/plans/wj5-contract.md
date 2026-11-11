# WJ-5 divergence diagnosis and local normal refinement

This work continues reviewed `6eb53fbde2437ee22e9597f4ef64ee8f6683a834`.
WJ-4 data, closed ledgers and its OBSTRUCTED horizontal qualification remain
immutable. The specification hashes are
`1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124` and
`0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c`;
worktree and original bytes agree. Main eq:integrals-linear retains the literal
`quad` typography erratum for `\quad`, with coefficient `(d+a*h)/a^2`.

## Fixed equations and boundary problem

All WJ-1 through WJ-4 model settings remain fixed: complete 3D tensor and cubic-N
surrogate, Original stress, DevStrainRawWorkV1, physical epsilon, Praw=-R:Graw,
independent Cd scalar diffusion; c1=1.8,c2=5/9,Ce1=1.44,Ce2=1.82,sigmak=1,
sigmaepsilon=1.3,Cd=572/6561,qe=81/(2 sqrt(143)),kappa=.381917863546733.
Molecular viscosity is omitted only in the outer high-Re model. Physical wall
nu=1.5e-5 remains positive. Linear full-strain inner stress and uncapped sources
are EP P_A-p_A,0 (chi=1,zm=H/2), or budget P_A-E_c (zK=zm=H,no physical chi).
K=k(H) remains an endpoint sample; budget inner k is transported and nonconstant.
Canonical u*=.5,H=1,top=2,h=.1,K0=.8466950351321795.

The domain remains [0,1]x[0,2], with reference differential-field Dirichlet
traces on both laterals and top, no-slip ground and zero modeled ground k flux.
Q3 velocity/scalars and regional continuous Q1 pressure, one overall gauge,
pi=p/rho+2k/3, shared U/k, separate pressure/gradient traces and one dependent
outer epsilon value trace remain unchanged. There is no inner epsilon PDE or
extra flux. IMEX1/PimpleSingleSweep, eta=0, zero prescribed forces/gravity,
reuse/verification Off, one U/pressure/k/epsilon solve and one logical
WJNativeV1 MassProjection operation are retained. StopOnBound=On; qualified
baselines require inactive numerical correction. The 3000-unknown guard stays.
No stabilization, stress/source clipping/capping, roughness smoothing, changed
boundary conditions, new backend or extra governing correction is introduced.

## Passive diagnosis and uncertainty

Retained WJ-4 startup and step32 states at Tx=.001024 use their actual producing
reports and accepted callback packets. A new observer lineage exports Gxx,Gzz,
divergence, both normal and horizontal face traces with explicit geometry,
actual B*U and lambda*m_p continuity rows, net flux and regional pressure mass.
Offline pressure-space L2 projection and element Legendre content diagnose the
Q3 divergence component invisible to continuous Q1 tests. They never correct U.
Actual mechanical pressure gradient and conservative tangential convection at H
form F_x; signed near-h storage/advection/stress/pressure/roughness loads use
unchanged production quadrature. Exact reference F_x is zero; viscosity is
smooth near h. No cause is inferred from weak/strong residual coexistence alone.

Constraint-aware uncertainty uses actual observation rows and q=q_D+P*q_F.
E_boundary*P is evaluated, including the native dependent epsilon trace action
on free endpoint K. Prescribed-data, lift, basis/geometry evaluation and
accumulation reserves remain, with the declared 100x arithmetic factor.
Pressure is unconstrained and retains a full-domain gauge-aligned estimate;
no regional gauge or fixed-velocity exemption applies to pressure. The original
WJ-4 boundary UNRESOLVED result is retained as a separate estimator outcome.
The dependent epsilon row propagates free k and free epsilon estimates through
their respective rows, with their physical units. Early/middle/end arithmetic
provides a conservative coefficient envelope; it is applied to every actual
startup/cadence observation map. Both compared trajectories retain the100x
reserve, including prescribed-data/lift and accumulation errors. A fixed value
trace does not fix its normal derivative: the tighter velocity-value estimate
cannot exempt strong divergence from its unchanged accuracy gate.

Producing/measurement source paths are frozen once used. A metadata or
measurement repair owns a separate checkout and versioned output lineage;
measurement-only qualification may reuse immutable native rows/packets without
another governing attempt. Development receipts remain visible. The native
state reader retains signed binary64 zero, while intent metadata retains JSON
integer types for canonical fingerprint verification.

## Frozen mesh/run design and stopping

Explicit validated z knots preserve defaults and bind exact x/z vectors,
comparison role, physical settings, dt, checkpoint and target in immutable intent.
M1 adds .085,.095,.125,.175 to all L1 knots. M2 adds to M1 .0825,.0875,.0925,
.0975,.1125,.1375,.1625,.1875. No knot is removed or physical interface moved.
Expected Nx1/Nx2 sizes: M1 1409/2447, M2 1713/2975.
Physical observations split actual compared meshes and physical/reporting
breakpoints; duplicate face owners have zero integration weight.

Per policy: one explicitly adverse L1/Nx2 dt=1.6e-5 64-step control to Tx;
M1/Nx2 dt=3.2e-5 32-step screen, M2 only if either M1 policy fails/unresolved.
Select one common normal vector using the worse outcome. If neither qualifies,
stop further evolution. Otherwise matched selected Nx1/Nx2 reach frozen
T3=.004096 at dt=3.2e-5 (128 total steps, reusing qualified screen CP32), and
selected Nx2 dt=1.6e-5 fresh CP0 reaches 256 steps. Accuracy cadence is32/64,
including startup and every common time; early/middle/end arithmetic is passive.
If a selected trajectory later fails, its target stays frozen and continuation
stops; no second campaign is authorized. A late32 segment is checked with an
extra16+16, immutable readback and one deliberate post-solve rejection.

All WJ-4 amplitudes/budgets and the 10% sensitivity criterion are unchanged.
Every required field/scope/time and independent evolved/weak-load check enters
qualification; missing/nonfinite/stale evidence cannot pass. The adverse
control is diagnostic, not a qualified baseline. Stationarity is not a target;
WJ-4's NOT_REACHED conclusion remains.

A NEW persistent WJ-5 ledger caps each policy at800 scheduled/conservative
control units and10800 numerical-study seconds, including failures, observations,
SVD, I/O, regressions and transactions. Unknown reservations remain charged.
Worst core is608 attempts plus33 transactions before controls. Costs are measured
in screening; remaining estimates are frozen before the selected full pair.
Linux is authoritative with isolated source-matched consumers, numerical threads1,
-ffp-contract=off -fno-fast-math, current capacity checks and <=72 total workers.
