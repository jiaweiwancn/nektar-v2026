# H26-A: prescribed-data native implicit flow contract

Starting reviewed/remote head: `10e3b1a60709314e1b93093fe03f6c8d77116bea`.
The branch owner is `.worktrees/wj-earsm-awf`; its untracked H25 pilot and the
main checkout's literature/work are preserved. H25/H25R1 accounts stay CLOSED.
H26 is a separate opt-in kernel and advancing executable, using native
`MultiRegions::ContField`, native geometry, quadrature, basis transformations,
physical derivatives and continuous local/global maps. It constructs sparse
mixed blocks, not the legacy dense regional pressure/mass machinery. The
kernel is exposed for subsequent EquationSystem/scalar integration; invoking
the dedicated H26 executable is the explicit opt-in. Legacy solvers stay intact.

## Equations, signs, time levels and spaces

Physical Cartesian G_ij=dU_i/dx_j; the two-dimensional computational x,y map
to physical x,z for the full three-dimensional WJ kernel. Retain c1=1.8,
c2=5/9, the physical N root, DevStrainRawWorkV1 and P=-R:G. Call the historical
pure WJ API with outerNu=0, then add physical molecular stress separately.
H26 retains positive nu in BOTH regions. Future scalar diffusivity uses the
independent Cd convention. No old high-Re log-law is declared an exact finite-nu
equilibrium and no balance forcing is supplied for its preservation.

M(U)=G+G^T, S=M/2, Pi=p/rho+2k/3; gravity potential is zero in H26 fixtures.
V=nu M-R+(2k/3)I. At a prescribed/current iterate *, choose positive
mu=nu+nu_t (the geometry-blended positive auxiliary/WJ coefficient in the wall),
Trem=V(*)-mu M(*). Verify V(*)=mu M(*)+Trem directly. Freeze ONLY this affine
stress in each linear problem. The nonlinear-flow fixture refreshes it at every
current iterate with prescribed positive k/E. There is no transported scalar.

For v with homogeneous essential velocity trace and q in continuous Q2:

    alpha (v,U) + (grad v, mu M(U)) - (div v,Pi) - (v,C_AWF Pi)
      = (v, f - EX2 convection + BDF history) - (grad v,Trem) + <v,t_total>.
    (q,div U) = (q,g).

Thus A=alpha M_u+K_sym, K_sym(v,u)=integral 2mu S(v):S(u), G=-D^T on the
closed test space, G_eff=G-C_AWF. On a natural boundary t_total=(V-Pi I)n,
which equals physical (-p/rho I+nu M-R)n. Do not add an old pressure jump.
All four velocity-component strain blocks are retained, including transpose
strain and off-constraint grad(div U). Full pressure coupling can be nonsymmetric.

Velocity Q3 (4 modes/direction), pressure global continuous Q2 (3 modes), common
GLL quadrature: 12x12 for tiny/curved verification, 8x8 for the directly matched
700-cell matrix (and higher-order independent quadrature convergence checks).
These are degrees, not names inferred from a pressure reduction convention.
The 700-cell physical mesh comes from successful H25 S0, unchanged x/height
knots, d0=.0875,H=.15,h=.10854503542535461,nu=6.423173803526447e-5.
The dimensional matching height remains 9.5 mm, Hw=6 mm, d0=3.5 mm. Small
verification boxes are declared manufactured fixtures, not recalibration.

Impose essential lifts in coefficient space before eliminating rows/columns.
Increment vectors have ZERO prescribed trace. Native GlobalSolve is never
called on an acceleration with inherited inhomogeneous inlet values. Native
local/global signs and topology are retained; no dense global Riesz inverse is
built. Pressure uses a mean multiplier only when all velocity boundaries are
essential (one genuine constant null mode); natural total traction fixes pressure
and uses no gauge. Independently test extra null modes in the tiny saddle matrix.

Constant-step BDF2: alpha=3/(2dt), history=(2U_n-.5U_nm1)/dt,
EX2 convection=2N(U_n)-N(U_nm1). Start using one genuine BE/EX1 step from the
declared initial fields, then retain the two accepted states/convection histories.
BE startup has O(dt^2) solution error; a targeted exact-history temporal fixture
can isolate BDF2 truncation independently. All viscosity, prescribed stress,
forcing, wall data and essential lifts are evaluated at NEW time. True restart
serializes both states, convection histories, pressure, time/step/dt, coefficient
rule and mesh/space identity; legacy histories are not imported as H26 data.

## Actual AWF pressure action and geometry

C_AWF Pi is +1_(r<h) E_H[P_H grad Pi|_H] on the RHS of strong momentum;
therefore its weak action is SUBTRACTED from G on the LHS. The known mechanical
pressure sample is grad Pi -2grad k/3 (minus grad potential if later retained),
so its known -2grad k/3 contribution belongs to prescribed body data, not G.
Use the same H sample and tangential projector as the existing affine action
in RANSWallEvolutionActions.cpp. It is not an approximate pressure Poisson law.

Flat geometry uses r=z-d0 and vertical extension. The curved manufactured graph
uses x=a,z=r+z_ground(a), a polynomial native curved coordinate map, and extends
vertically at the same a. Here r is the declared graph coordinate, not Euclidean
distance; physical projector n=(-z_ground',0,1)/sqrt(1+z_ground'^2), physical
gradients/Jacobians and volume/surface measures are used. This is a specified
general-flow extension, not a claim of hill calibration or exact normal-distance
coordinates. Split the roughness integral at r=h; do not sample an unsplit
indicator on ordinary GLL nodes. Both native trial and independent oracle use
the actual pressure trace polynomial on H, with its physical tangential gradient.

## Prescribed wall tensor adapter

Use P=I-nn^T, b_A=-nu_A P M n,
Rh=P RW P+(n.RW.n)nn^T+b_A n^T+n b_A^T,
Rwall=(1-chi)Rh+chi RW. Use geometry-only chi=0 for r/H<=.4, 1 for r/H>=.75,
and quintic smoothstep 6s^5-15s^4+10s^3 between. This is C2 and identically full
WJ in the open collar (.75H,H) and above. The same limiting prescribed k,E_w/
epsilon and physical gradient are supplied on both sides; no C1 constraint is
imposed on native C0 velocity gradients. Test full physical weak traction and
convergence separately. Complete blended tensor enters weak load, so derivatives
of chi and geometry are retained by integration by parts, not omitted.

Planar oracle: physical root N>=c'=1.8 solves
N^3-c'N^2-.35q^2 N-c'q^2=0, q=kG/Ew,
C=3N/[5(N^2+q^2)], d=Cq^2/N.
Rh_xr=-nu_A G; Rh_xx/k=2/3+d, Rh_yy/k=2/3, Rh_rr/k=2/3-d.
Collar shear is -[(1-chi)nu_A+chi C k^2/Ew]G. Without full-WJ collar, equality
of endpoint shear would additionally require nu_A(H)=C(qH)K^2/epsilonH.
Trace 2k, rotational covariance, stress work, realizability and collar recovery
are gates; no clipping. Revised cutoff/reference budget and scalar admissibility
are H27 work. The old monotonicity argument is not reused for cutoff-dependent Ew.

## Linear algebra and coefficient refresh

Retain native element basis/geometry/mapping, sparse mass/D/G/C and graph. Form
current K_sym from native physical derivatives. Initially use a sparse mixed
factorization as a retained right preconditioner, and restarted GMRES on the
CURRENT mixed operator. A previous factorization is never the physical operator.
Recompute true momentum/continuity residuals after solves. Refresh preconditioner
when coefficient relative range changes >25%, iteration count exceeds 30 or
twice its initial warmed value; record every refresh. Cap Krylov at 400 actions.

For declared affine coefficient space mu=nu+a(t)mu_t(x), precompute the two
native K bilinear forms and combine at current a(t), exactly by linearity. This
is current implicit viscosity, not an explicit difference term. Independently
compare against fresh quadrature assembly after a material coefficient change.
General/nonlinear mu is reassembled from current quadrature data. Profile both
assembly choices; record native variable-coefficient cache/lifetime effects rather
than relying on undocumented cache safety. Tiny direct matrices are oracles only.

## Frozen validation scales and fixtures

Dimensionless physical/reference scales: Lx and Lz from each fixture, U0=1,
P0=.01, k0=.008809015145515192, E0=.002315524054799254,
mu0=nu+.3819178635467331*.051*Lz/2. Use volume-normalized L2 velocity/pressure
errors (remove only the declared null mean), H1 velocity seminorm and native
mass-dual momentum/pressure-continuity residuals, scaled by reference inertia/
diffusion and U0/L. Report null RHS/action separately, with roundoff bound based
on actual operator and lift scales; no dimensionful max(1,||rhs||).

Gates: independent operator/traction action relative error<=1e-10;
independent normalized momentum AND continuity residual<=1e-9;
tiny direct-saddle solution discrepancy<=1e-8. Tolerances are not loosened after
results. Pressure gauge/nullspace rank is tested independently by SVD. Algebraic
tolerance target 1e-11, with one 1e-12 tightening check subordinate to spatial error.

Flat and curved tiny tests include nonconstant mu, both components, cross strain,
nonzero pressure, nonzero essential lift and a natural total-traction boundary;
an active-AWF case needs nonzero pressure correction. Verify dissipation/work,
material mu change with retained preconditioner, adapter and nonlinear frozen-
scalar refresh. Deliberately missing transpose strain, flipped C sign and stale
mu actions must fail their discrepancy tests.

Independent Python oracle evaluates Modified_A polynomials/Jacobi derivatives,
analytic graph Jacobian and separately refined Gauss quadrature using exported
native topology/signs only. It assembles a tiny direct saddle problem. Continuous
MMS forcing is generated symbolically from physical derivatives, not A_h u_h+G_h p_h.
Use polynomial streamfunction for temporal isolation (Q3 velocity/Q2 pressure),
and a trigonometric smooth streamfunction for spatial refinement. At least three
dt levels and three h levels; expected affine smooth orders U L2=4, H1=3,
pressure L2=3; curved geometry order sufficiently high to resolve mapping error.
Require observed BDF2 velocity order>=1.8; declare pressure temporal target>=1.8
when spatial/splitting errors are resolved. Report weak divergence AND strong
divergence; do not demand pointwise incompressibility of a weak FE solution.

## Matched controlled engineering screen

Cases Snu/Smu are actual stock no-RANS native solvers at nu/mu_ref; Nnu/Nmu
disable extra stress and AWF. Nspace uses positive spatial viscosity/stress;
Ntime modulates the turbulent part smoothly by10% and has active C_AWF.
Use representative mu_t from retained admissible reference data (approximately
.3819178635467331*.051*r outside, auxiliary positive extension inside); pin exact
range and mu_ref from generated inputs before benchmark execution.

All share materialized700 mesh, Q3/globalQ2, GLL8, one rank, six library thread
limits1, IMEX/BDF2 order2 with BE startup, dt .005, same physical intervals,
warm-up and boundaries where stock supports them. Native full-strain/coupled
pressure versus stock component/split pressure is an explicit remaining difference.
Natural total traction is qualified separately when stock cannot match it.
Manufactured forcing is verification-only. No old log-law preservation forcing.

Preregister 20 warm-up advances then400 accepted advances per block, three
blocks per case, interleaved by repetition. Each warmed block must exceed3 s;
extend ONLY if timer dominance or variation leaves the decision ambiguous.
Record mean/median/spread, iterations/corrections, active C_AWF Pi, coefficient
changes, rebuilds and every ordinary source/BC/residual/logging cost. Require
time-case pressure feedback and correction norms to exceed1e-6 of their declared
reference scales in warmed blocks, not merely a configuration flag.

Freeze recovery every200 accepted steps, compact scientific summary every20,
and final output. Checkpoint serialization and restart read/restore are separately
timed; bound recovery retention to latest2 plus protected test/terminal states.
Disjoint categories: coefficient/constitutive/source updates; equation assembly/
application/solution; preconditioner setup; true residual/diagnostics; I/O; remaining
controller overhead. Attribution>=95%; no inclusive timer double-counting.

Tstar=505.02466521278404, dtstar=.005, Nstar=101005.
Cscreen=1.25[Cfixed+101005*c_flow+Crecurring]<=43200 s.
c_flow is the larger largest warmed block mean of Nspace and Ntime. Cfixed is
measured startup+one actual restart read/restore+final output. Count EVERY
scheduled checkpoint/summary not already inside ordinary means exactly once.
This is a half-day kernel screen, not a production/stability/hill claim.

## Resource and checkpoint plan

New ext4 root `/home/jay/build/wj-h26-native-flow-20261117`, no old capsule or
ledger modification. Initial Tailscale inventory:37 workers, ~72GiB free, ext4;
USB timed out. General Linux checkout f4b70bae with unrelated developer-guide
change is preserved. Initially serial science, build at most8 jobs, worker
admission/monitoring includes helpers and all unrelated jobs, hard ceiling72.
RSS8GiB science/16GiB build;8GiB source/build and8GiB science;15GiB free reserve.

New aggregate allowance21600 operation seconds;1800 protected for evidence and
publication until their operations are charged. Initial control inventory/setup
is conservatively fully charged60 s; subsequent bounded operations measured,
failures/unknowns retained, overlaps separately charged. Initial pilots<=120 s,
build<=900 s, numerical campaigns adapt from pilot costs with per-run limits.
Before duplicate launch resolve any unknown outcome; never stop unrelated work.

A publishes contract, independent fixtures and namespace/admission/account.
B publishes compiled advancing core and actual numerical/preliminary profile
evidence, including failures and subsequent repairs. C publishes final selected
candidate/barrier, all identities/costs/accounts and H27 contract. All are pushed
and exact remote commits verified. H26_NUMERICS, H26_COST and readiness are
separate; FULL_COUPLED_MODEL_AND_PRODUCTION remains NOT_QUALIFIED.

## Prospective implementation refinement before the controlled matrix

The v3 complete-step pilot measured 0.703085529 s/step (four warmed steps):
source/projection 0.212 s, equation application/solution 0.486 s, and diagnostics
0.0054 s. It is a rejected slow prototype, not the engineering decision. The
v4 pilot measures 0.0704072744 s/step (20 warmed steps); the controlled three
blocks per case remain required. No acceptance target or equation changed.

Retain exact affine sparse coefficient/operator bases and analytic spatial
source factors. For U=e^-t U0, Pi=e^-2t Pi0 and mu=mu0+a(t)mu1, continuous
manufactured f=e^-t f1+e^-2t f2+e^-t a(t)f3. Each spatial load is integrated
from continuous physical derivatives at startup, never generated by applying
the discrete saddle operator to a manufactured state. Current a(t), BDF2
history, EX2 convection, exact new-time boundary lifts, true residuals and
normal logging are included in every ordinary step. Startup includes caching.
Native current matrices are exact affine combinations, with a separate BE
mass coefficient at the first step. The physical mean-pressure operator is
unchanged; only its sparse preconditioner factors a pressure-pinned equivalent
and projects the inverse back to the mean gauge on every application.

The actual stock no-RANS comparator uses unchanged VCSWeakPressure and native
IMEX(2), with its usual pressure Poisson and component Helmholtz solves. A
prescribed-force callback makes the stock explicit source combination equal
the same new-time continuous body force as the candidate. For rho=exp(-r dt),
the callback at time n dt returns w_n=A rho^n+(rho-A)2^-n, A=rho²/(2rho-1),
for each exponential factor. Thus w_0=rho gives genuine BE startup and
2w_n-w_(n-1)=rho^(n+1) gives the exact BDF2 new-time source. This changes no
stock solution operator and introduces no RANS data. Source cache construction,
callback work, BCs, field transforms and current diagnostics are timed. Prove
this scalar identity independently and verify constant-limit continuum errors.

All six cases use the same polynomial continuous velocity/pressure, BE first
step, BDF2/EX2 thereafter, dt=.005, 420 accepted steps over [0,2.1], warmup20,
400 retained steps, common XML geometry/degree/quadrature/boundary traces,
one MPI rank and six numerical-thread controls1. The stock pressure splitting
and component strain versus the candidate's converged full-strain mixed
coupling remain declared differences; their continuum constant-mu limits
agree on divergence-free fields. Full total-traction boundaries are qualified
on separate manufactured fixtures. Stock comparison-state outputs have the
common cadence/shape, but are not claimed to restore the stock GLM/HOPBC
histories. Only the candidate's actual complete-history restart is qualified.

The nonlinear adapter cache preserves element/quadrature identity. Native C0
gradients from coincident points in adjacent elements must remain distinct.
The v3 coordinate-keyed prototype failed the complete-stress weak identity;
v4 fixes that bug and passes it at 3.31e-17 without imposing a C1 constraint.
