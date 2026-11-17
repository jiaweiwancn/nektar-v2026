# WJ-H27 — Native finite-nu scalar and wall coupling handoff

H26 supplies a numerically verified serial native Q3/global-continuous-Q2 flow
kernel, current full symmetric-strain diffusion, the actual nonlocal AWF pressure
action, a prescribed wall/WJ tensor adapter, and complete flow-history restart.
Its selected source is `8aeb00a2f37437cd17e043a9ba6f86a6f208e1d6` on
`codex/wj-earsm-awf`. The final H26 report and evidence index give the executed
performance screen. Full coupled model and production remain NOT_QUALIFIED.

## Coupled mathematical contract to establish

Advance velocity/modified pressure together with transported inner k, outer k
and outer epsilon. Inner E_w is dependent wall dissipation, not a separately
transported epsilon field. Compute current K and the outer epsilon trace at H,
and derive E_w and the revised cutoff/reference budget from the current state.
Retain physical molecular viscosity in the new outer momentum and revised
finite-nu scalar equations. Keep the pure historical high-Re WJ kernel at
`outerNu=0`; add molecular terms in the new PDE operators explicitly.

Use full WJ with c1=1.8, c2=5/9, the physical N root and
`DevStrainRawWorkV1`. Evaluate P_turb=-R:grad(U) from the complete current tensor
and raw physical gradient. The general three-dimensional explicit formula uses
the documented cubic-root surrogate; it is not an exact solution of a general
three-dimensional implicit ARSM. Do not replace it with the planar simple-shear
formula in a hill flow.

At every current nonlinear iterate use

```
V = nu*(grad U + grad U^T) - R + (2*k/3)*I
mu > 0
Trem = V - mu*(grad U + grad U^T)
G_eff = G - C_AWF
```

The current implicit operator must contain the current positive mu and all
transpose-strain/variable-coefficient terms. A retained old preconditioner may
approximate that system; it must not become the physical operator. If Pi
contains 2k/3, sample mechanical pressure with the known `-2*C_AWF*k/3` term,
using a declared trace/interpolation space for transported k. Declare gravity
and all scalar diffusion/source signs and time levels before implementation.

Retain the proposed adapter

```
P = I - n*n^T
b_A = -nu_A*P*(grad U + grad U^T)*n
R_h = P*R_W*P + (n^T*R_W*n)*n*n^T + b_A*n^T + n*b_A^T
R_wall = (1-chi)*R_h + chi*R_W
```

Use a geometry-only smooth chi with an open full-WJ collar below H and common
limiting physical k, E_w/epsilon and gradient inputs. H26's quintic transition
on r/H in [.4,.75] is one tested choice. Identical constitutive formulas do not
make C0 velocity gradients continuous: retain weak physical total-traction and
refinement tests. Keep complete blended-tensor weak loads, geometry terms and
all chi derivatives if strong differentiation is introduced. Do not impose an
old pressure-jump law or silently clip stresses to repair realizability.

Re-derive reference production, finite-nu velocity anchors, cutoff admissibility
and existence/uniqueness. The new reference production can depend on the cutoff
through E_w. The old monotonicity argument and old cutoff/reference cache are
not a proof for the revised model. Declare failed/admissible branches and their
continuation policy. Treat the finite-nu log-law inlet as boundary data, not as
an exact steady state maintained by hidden balance forcing.

## Native implementation entry points and limits

`EquationSystems/WJNativeFlow.{h,cpp}` provides signed native continuous maps,
native quadrature/physical derivatives, sparse mass/Riesz operators, all strain
blocks, physical Dirichlet lifts, a nonsymmetric mixed GMRES solve, and true
mass/dual residuals. `StrainSamples` and `StressSamples` accept elementwise
quadrature data. Preserve element identity at coincident points: C0 gradients
on adjacent elements generally differ. The failed coordinate-keyed prototype
and its corrected complete-stress weak identity are retained in H26 evidence.

The sparse pressure-pinned factorization is only a preconditioner. Its inverse
is projected back to the physical mean gauge; the current target operator still
contains the mean constraint. Natural total traction removes the pressure
nullspace in the tested fixtures. Check extra null modes whenever spaces or
boundary laws change. Mass/Riesz solves for acceleration must use homogeneous
increment spaces, independently of nonzero velocity inlet data.

The H26 affine operator/source caches are exact for the prescribed test data.
General evolving WJ/scalar fields need fresh current quadrature coefficients
and reassembled or matrix-free current actions. Retain topology, geometry and
useful preconditioners, measure coefficient/iteration changes, and validate any
refresh/cache policy against uncached current actions. The default policy
refreshes for a supplied relative coefficient change above25%, more than30
linear iterations, or more than twice the initial iteration count. Retune only
from measured coupled data while preserving final residual targets.

The dedicated H26 driver is opt-in and serial. Integrate an explicitly named
native coupled EquationSystem/factory path without changing historical solver
defaults. Qualify distributed signed assembly, nonlocal H sampling, pressure
gauge handling and restart before using multiple MPI ranks. H26's structured
graph fixtures use r=z-ground(x) and vertical E_H extension, with physical
tangential projection and measures. A general hill/unstructured mapping needs
an explicit geometry/point-location/communication contract and new tests.

## Required H27 qualification

Keep the H26 operator/direct, known-k pressure sign, tensor, work/traction,
current-coefficient, finite-nu continuous-MMS, BDF2, tolerance and restart gates.
The nonlinear flow fixture must have zero continuity forcing: the predecessor
fixture's interpolated divergence datum did not close that gate. H26's corrected
native mixed projection and final zero-continuity nonlinear solve pass it.

Add independently differentiated coupled continuous manufactured solutions for
velocity, pressure and positive evolving scalars, with spatial/time refinements
and nonzero pressure/AWF/known-k actions. Verify scalar conservation, interface
fluxes, stress work and revised budget residuals. Check constitutive/scalar
admissibility and report failures without unlabelled clipping. Test the complete
coupled current residual after refreshing U, k, epsilon, K, E_w, cutoff and
stress data; a frozen-scalar flow residual does not bound coupled convergence.

Extend restart to every scalar time level, explicit convection/source history,
wall/cutoff state and any branch/hysteresis state. Pin all scientific parameters,
geometry, time-integrator data and implementation contract. Compare resumed and
uninterrupted coupled trajectories at the same completed physical time, and
reject mismatched scientific inputs. Historical H25/H25R1 histories remain
diagnostic artifacts, not qualified revised-model initial states.

Measure coupled complete-step work, nonlinear iterations, setup/refreshes,
admissibility/diagnostic work and actual output/recovery costs separately. H26's
six-case frozen/prescribed screen and small nonlinear test do not forecast
transported-scalar coupling cost. Use a new finite resource allowance; the H26
ledger is closed after its final publication and must remain immutable.

## H28 boundary

H28 qualifies the full model's hill mesh/order/timestep, nonlinear robustness,
physical convergence, monitoring, restart/storage recipe and one-day production
forecast. The historical T_star=505.02466521278404 and dt_star=.005 are H26
engineering yardsticks, not qualified revised-model production parameters.
Keep the shared-host total at most72 actual workers, explicitly budget rank and
library threads, and store heavy data in an isolated approved Linux namespace.
