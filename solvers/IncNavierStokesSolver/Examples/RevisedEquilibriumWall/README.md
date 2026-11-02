# Revised CLS/EP-AWF native wall stage and mapped evolution


`SolverType=RANSEquilibriumWallEvolution` selects the M4 time-dependent
multiregion solver. It retains the original CLS constitutive law and the
revised-note physical EP-AWF source. The two selectors share the same declared
flat physics; the original M3 stage remains available for its independent
spatial checks.

| Selector | Operation | Qualified mesh and execution |
| --- | --- | --- |
| `RANSEquilibriumWallSystem` | One supplied nonlinear stage; `NumSteps=0` | Serial 2D affine quadrilaterals, M3 |
| `RANSEquilibriumWallEvolution` | Native IMEX1/2 with native startup and full-history restart | Affine charts or analytic mapped quads/hexahedra; serial and MPI 2/4, M4/M5 |

M4 uses globally numbered native CG coefficients replicated on each rank.
Physical cells own their volume/source quadrature through `geometry ID % ranks`;
loads are reduced, and rank zero solves the dense coupled system and broadcasts
the update. This bounded verification backend is limited to 3,000 coupled
unknowns for nonlinear solves. Larger meshes permit reference-only audits.
It establishes partition agreement for the tested flat cases; it is not a
scalable distributed sparse implementation. Supply the original unpartitioned
XML mesh when launching MPI. Generated evolution sessions use the parallel-compatible
world setting `GlobalSysSoln=IterativeStaticCond`; rank-local native CG sessions
use `DirectFull`, and the coupled wall solve uses its own dense backend.

## M4 reference and native time integration

The [unanchored reference fixture](../../Tests/RANS/EquilibriumWallReference.xml)
uses chi=1, eta=0, p=3/p=1, H=2, Z=30, length 300, roughness 0.1,
r_K=r_m=1, nu=1.5e-5 and u_star=0.4. It audits the initial state with
`NumSteps=0`. The [evolution smoke fixture](../../Tests/RANS/EquilibriumWallEvolution.xml)
uses the identical projected reference and eta=1 for three actual IMEX2 steps.
Both use two horizontal cells so the free-load constant tests and energy work
are nontrivial. A 2D mesh uses physical (x,z); a 3D mesh uses physical (x,y,z)
and variables `u,v,w,k,eps_tilde,p`.

The finite-dimensional target is a constrained physical L2 projection of the
composite AWF/logarithmic velocity and outer A/z dissipation. It represents
constant k and horizontal invariance exactly. The target also preserves no-slip,
the matching velocity, its one-sided reference gradients, the dependent epsilon
trace, and compatible top values/gradients. These are constraints on the stored
native polynomial, not replacements for its evaluated derivatives. Inner and
outer pressures have the offset -11*K0/75 and one common volume-mean gauge.
The validator independently reconstructs the constrained L2 minimizer, gradients,
full matching traction, zero k flux, top epsilon flux and gauge.

`RANSWallChi` and `RANSWallEta` are independent. Chi=1 selects physical EP-AWF
redistribution. Eta=0 evaluates that model without an extra discrete correction.
Its polynomial target generally has nonzero differential residuals; the
validator records those defects first and measures h/p refinement. Eta=1 is
restricted to isolated flat neutral verification. It subtracts the one stored
free momentum/k/outer-epsilon load, after checking the unchanged divergence,
gauge, essential and dependent-field constraints. It never corrects those
constraints or resets free state coefficients. The JSON reports the complete
unmasked reference loads, removed rows, fixed free defects, their constant-test
integrals and energy work. This subtraction changes the discrete equations.

Declare the native scheme explicitly:

```xml
<TIMEINTEGRATIONSCHEME>
  <METHOD>IMEX</METHOD>
  <ORDER>2</ORDER>
</TIMEINTEGRATIONSCHEME>
```

Only orders 1/2 with an empty variant are supported. IMEX2 uses the unchanged
native first-order startup. Every state-dependent transport, stress, reaction,
sampler and pressure/constraint coupling is implicit. Prescribed, state-independent
body accelerations are explicit, including their prescribed part of the
roughness force. The full residual and native explicit derivative share one
fixed physical body-load quadrature, independent of the current sampled K;
state-dependent roughness forcing retains the adaptive wall quadrature. Additional manufactured scalar sources are evaluated at the
implicit stage time. The differential state is flattened only for the native
time-integration interface; each field retains its own physical mass matrix,
space, volume and constraints. Pressure and the common gauge remain algebraic
stage unknowns.

With native combined history h and stage coefficient lambda, the free rows solve
`M*(Q_stage-h)/lambda + R_implicit(Q_stage,t_stage) - eta*d = 0`, together
with the original mass, gauge, boundary and dependent-field constraints. The
native explicit history supplies the body-load contribution to h. At the
zero-body flat target, consistent value histories equal Q0, derivative histories
vanish and h=Q0; the stored d cancels the same stage's spatial differential
balance. The pressure and boundary path is this coupled variational problem,
with no separate velocity-correction pressure solve. Acceptance uses the actual
full residual; there is no final reference projection. The fixtures explicitly
use the scaled physical stage tolerance 2e-12, stricter than the 1e-9 reference
residual/state/constraint gate.

| M4 selection or parameter | Meaning |
| --- | --- |
| `RANSWallReference` | `Neutral` (default) constructs the constrained target; `Functions` uses analytic initial/boundary/manufactured data |
| `RANSWallEta` | Global 0 (default), or 1 only with `Neutral`; no spatial eta |
| `RANSWallFrictionVelocity` | Positive reference friction velocity, default 0.4 |
| `RANSWallDirectionX`, `RANSWallDirectionY` | Fixed unit horizontal direction; spanwise flow requires a 3D mesh |
| `Time`, `TimeStep`, `NumSteps` | Initial reference epoch, fixed positive dt, and number of additional native steps |
| `RANSWallCondense` | Eliminate only free inner velocity; retain inner k and every pressure/gauge/history mode |
| `RANSWallRestart` | Full checkpoint directory, or `None` |
| `RANSWallCheckpoint` | Output directory base, or `None`; unique generations prevent overwriting |
| `IO_CheckSteps` | Optional checkpoint interval; the final state is also checkpointed |

`Functions` mode requires inner/outer initial velocity/k/pressure, outer initial
epsilon, external velocity/k boundary functions, and outer epsilon boundary
functions. Epsilon at H is still reconstructed from current sampled K. Source
function names match the M3 table below. `StageHistory` functions are not used:
the native integrator owns all time histories. The generator
[cls_epawf_m4_cases.py](../../Tests/RANS/cls_epawf_m4_cases.py) supplies reproducible
2D/3D neutral and 2D evolving manufactured sessions.

Run a copied M4 fixture in a fresh Linux run directory using the environment
and isolated executable shown below. The complete qualification is registered
as `RANSCLSEPAWFM4Oracle`:

```bash
ctest --test-dir "$CLS_EP_RUN/build" --parallel 1 --output-on-failure \
  -R 'RANSEquilibriumTimeNative|RANSEquilibriumEvolutionNative|RANSEquilibriumWallEvolutionSmoke|RANSCLSEPAWFM4Oracle'
```

The oracle runs its MPI cases sequentially with 1/2/4 ranks and one numerical
thread per rank. It checks all three chi values in independent 2D/3D residual
oracles, h/p refinement, constrained target and anchor budgets, unanchored native
steps, a nonpolynomial prescribed-load partition control, several dt values,
one/repeated steps, full 3D horizontal rotation, nonstationary manufactured
orders, full/condensed agreement, MPI assembly/solve agreement, and full-history
restart across startup and partitions. The native time probe separately checks
the actual IMEX recurrence, every history slot, native startup and invalid
checkpoint rejection.

A checkpoint contains the coupled state, pressure/gauge, all native value and
explicit-derivative slots, their canonical time labels, fixed reference/defect,
eta, and a mesh/basis/physics/function signature. `State.bin` is checksummed and
published atomically in a new directory. A restart session keeps the original
`Time` epoch and the same physical/functions/scheme/dt data; the actual
continuation time and step number come from the payload. `NumSteps` counts
additional steps. Source, shape, scheme, dt, corruption and incomplete-publication
mismatches are rejected. Full-history continuation is compared with an
uninterrupted run in the same gauge.

The `.evolution.json` file records actual native output/history, stage and step
norms, field maps, cell ownership, anchor budgets and restart provenance. Native
field outputs retain their separate physical regions, including inner velocity
and pi_w. Large outputs and checkpoints belong in the Linux run directory.
The [M4 report](../../../../docs/plans/cls-epawf-revised-m4.md) records the tested
revision, commands, numerical results and retained development attempts.

## M3 single-stage reference

`SolverType=RANSEquilibriumWallSystem` selects the revised-note wall element.
The implementation uses native continuous spectral elements for velocity and
k across the full domain, outer-only epsilon, and separate inner/outer reduced
pressure spaces. It applies one global volume-mean pressure gauge. Mesh
coordinates `(x,y)` represent physical `(x,z)` and variables `(u,v)` represent
`(U,W)`. The outer CLS stress retains all three modeled normal components.

This M3 backend solves one supplied nonlinear stage in serial on affine,
axis-aligned, canonically ordered quadrilaterals. The mesh must have a complete
element boundary at the constant matching height H. Polynomial bases are
Modified_A, with equal degrees in both directions and each field's degree at
least two. The supplied fixtures use velocity/scalar degree p and pressure
degree p−2, with p at least four.
The dense coupled backend accepts at most 3,000 unknowns. Unsupported meshes
and MPI partitions are rejected by this M3 selector. The M4 selector below adds
native IMEX stepping, restart and full three-dimensional flat fields. Curved
geometry remains M5 work.

The [stationary manufactured fixture](../../Tests/RANS/EquilibriumWallFlat.xml)
has two elements, p=4, H=2, top height 4, roughness 0.1, r_K=r_m=1 and
nu=1.5e-5. Explicit manufactured sources maintain positive constant turbulence
at zero velocity. The initial velocity is perturbed, so the fixture exercises
an actual nonlinear solve. It is a verification case, not an ABL simulation.

Run on the Linux server after building the complete external solver:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
CLS_EP_RUN=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised
CLS_EP_SOURCE=$CLS_EP_RUN/source
export TMPDIR=$CLS_EP_RUN/tmp
export LD_LIBRARY_PATH="$CLS_EP_RUN/build/geometry-lib:$LD_LIBRARY_PATH"
mkdir -p "$CLS_EP_RUN/manual-flat" "$TMPDIR"
cp "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml" \
  "$CLS_EP_RUN/manual-flat/case.xml"
cd "$CLS_EP_RUN/manual-flat"
taskset -c 0 "$CLS_EP_RUN/build/IncNavierStokesSolver" case.xml
```

Use a fresh run directory. Keep generated fields, JSON, `.opt` files and logs
outside the source checkout. See the [implementation plan](../../../../docs/plans/cls-epawf-revised-implementation-plan.md)
for the isolated build commands and the 72-worker server limit; its initial
build allocation is 16 workers, with tests run sequentially.

The stage equation contains `alpha * (state - StageHistory)` for both inner
and outer velocities, full-domain k and outer epsilon. `RANSWallStageAlpha`
supplies alpha; `Time` supplies the evaluation time. `NumSteps` must be zero.
The JSON explicitly records `native_time_map=false` because this operation
does not invoke a native time-integration scheme.

| Session selection or parameter | Meaning |
| --- | --- |
| `TurbulenceModel` | `CLSHighRe` (default) or `StandardKEpsilon` control |
| `RANSWallCondense` | `False` (default) or `True`; eliminate free inner velocity modes, retain every pressure mode |
| `RANSWallStageMode` | `Solve` (default) or `Residual`; the latter audits the projected input without solving |
| `RANSWallHeight`, `RANSWallRoughness` | Constant physical matching height and roughness |
| `RANSWallSampleHeight`, `RANSWallTailStart` | Independent physical r_K and r_m |
| `Kinvis` | Physical inner molecular viscosity; absent from the declared outer high-Re operator |
| `RANSWallChi` | Fixed correction weight, default 1 |
| `RANSWallStageAlpha` | Nonnegative stage storage coefficient, default 0 |
| `RANSWallTangentialOrder`, `RANSWallNormalOrder` | Physical Gauss orders, defaults 8 and 16 |
| `RANSWallNonlinearTolerance`, `RANSWallMaxIterations` | Default coefficient-residual tolerance 2e-12 and iteration limit 15; physical block norms are reported separately |

All session functions use analytic expressions. Required initial functions
are `InitialConditionsInner` for u,v,k,p and `InitialConditionsOuter` for
u,v,k,eps_tilde,p. Optional `StageHistoryInner/Outer` have the same variables;
when omitted, history is the initial state. Optional `ExactSolutionInner/Outer`
support field-error reporting.

`BoundaryValuesInner` specifies u,v,k on the external side faces;
`BoundaryValuesOuter` specifies u,v,k,eps_tilde on the external sides/top.
Ground velocity is strongly zero and ground k flux is naturally zero.
The shared velocity/k traces are native CG traces. Epsilon at H is computed
from the current physical sampled K, so side epsilon corner values must agree
with that dependent trace. No inner epsilon PDE or independent pressure
Dirichlet boundary is created.

Optional `BodyForceInner/Outer` specify u,v body accelerations. Optional
`ScalarSourceInner` (k) and `ScalarSourceOuter` (k,eps_tilde) are additional
physical manufactured source densities. The inner tangential roughness force
includes `(1-a_h)*F_A`; storage and residuals remain unmasked. F_A uses the
one-sided outer value of `d_x(U^2)+d_x(pi-2k/3)-b_x` at H.

The solver writes eight native field files with their own element sets:
outer u,v,eps_tilde,p; full k; and inner u_w,v_w,pi_w. Both pressure outputs
use `pi=p_mechanical/rho+2k/3`. The `.patch.json` file retains the complete
stage state/history, residual blocks, physical mass-dual norms, constraint defects,
sampled K and forcing. Inner histories remain present when condensation is on.
The independent validator also reports residuals, state errors and storage
errors with the physical scales frozen in the M0 contract.

To generate an evolving manufactured case, import `evolving_session` from
[cls_epawf_m3_cases.py](../../Tests/RANS/cls_epawf_m3_cases.py) and call it with
the fixture path, destination path, degree 6, time 0.3, mode `Solve` and the
desired `condensed` boolean. NumPy/SciPy support the independent oracle and
SymPy constructs its analytic forcing. The complete reproducible qualification is:

```bash
taskset -c 0 python3 \
  "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m3.py" \
  --build "$CLS_EP_RUN/build" --output "$CLS_EP_RUN/manual-m3-qualification"
ctest --test-dir "$CLS_EP_RUN/build" --parallel 1 --output-on-failure \
  -R 'RANSEquilibriumPatchNative|RANSEquilibriumWallStageSmoke|RANSCLSEPAWFM3Oracle'
```

The original CLS adverse-diffusion and inner-k positivity findings remain
failed physical gates. A converged flat stage does not qualify long-time
stability or terrain accuracy.

## M5 analytic normal-coordinate geometry

`RANSEquilibriumWallEvolution` also supports analytic mapped quads/hexahedra.
The XML mesh represents an affine computational chart `(q1,q2,rho)`; in 2D
its native `(x,y)` means `(q1,rho)`. Physical positions are
`X(q)+rho*H(q)/H0*n(q)` and the matching chart is `rho=H0`.
Velocity components remain Cartesian. All physical gradients, volume masses,
pressure means and oriented face fluxes use this map, including the outer tube.
This path supports the analytic maps below; arbitrary Cartesian curved XML
imports use a different geometry representation.

| Selection | Physical base wall |
| --- | --- |
| `RANSWallMap=Flat` (default) | Horizontal wall |
| `RANSWallMap=Plane` | Graph `z=SlopeX*q1+SlopeY*q2` |
| `RANSWallMap=Cylinder` | Outward circular cylinder of `Radius`, extruded in q2 |
| `RANSWallMap=Gaussian` | Smooth graph `Amplitude*exp(-(q1/LengthX)^2-(q2/LengthY)^2)`; omit q2 in 2D |
| `RANSWallMap=CosineHill` | C1 graph `Amplitude/2*(1+cos(PI*q1/LengthX))` inside the feet, zero outside |

Map parameter names have the `RANSWallMap` prefix: for example
`RANSWallMapRadius=12` and `RANSWallMapAmplitude=.025`.
`H(q)=RANSWallHeight+HeightSlopeX*q1+HeightSlopeY*q2+
HeightAmplitude*sin(HeightWaveX*q1)*cos(HeightWaveY*q2)` uses the same prefix
for its height coefficients. Roughness, sample height r_K and tail start r_m
are physical normal distances. Each must fit below the certified minimum H.
The complete inner and outer normal tube must pass the positive-stretch and
injectivity certificate. A cosine foot at q1=±LengthX requires an explicit
native element/patch face; an element crossing it is rejected.

`RANSWallFunctionCoordinates=Physical` is the default for analytic initial,
boundary, body and source functions. It passes actual physical `(x,z)` in 2D
and `(x,y,z)` in 3D. `Chart` passes `(q1,rho)` or `(q1,q2,rho)` and is useful
for the supplied strong manufactured cases. In both conventions, vector
function values specify Cartesian components, including the Cartesian
vertical component of a tangent force on an inclined surface.

Use `RANSWallChi=0` and `RANSWallEta=0` for curved or variable-H M5 cases.
Nonplanar chi and curved anchoring are rejected. The isolated constant-H
inclined-plane reference can use chi=eta=1 for rotation/fixed-point checks.
M5 supplies geometry and conservative coupling; it does not select a spatial
transition or implement M6 sensors.

Generate a small circular manufactured run on Linux, outside the source tree:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
mkdir -p /media/jay/新加卷/nektar-v2026-runs/m5-circle-example
cd /media/jay/新加卷/nektar-v2026-runs/m5-circle-example
python3 - <<'PY'
import sys
from pathlib import Path
source=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/source')
sys.path.insert(0,str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
from cls_epawf_m5_cases import mapped_session
mapped_session(source/'solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml',
               Path('case.xml'),family='circle',degree=5,swirl=True,
               steps=2,order=2,dt=.005)
PY
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/build/IncNavierStokesSolver case.xml
```

`family='variable'` generates a smooth variable-H graph; `surface3d` with
`dimension=3` exercises both principal curvatures and both height derivatives.
`feet` supplies the explicit C1 joins. These are verification data with strong
manufactured sources, rather than physical terrain-flow validation cases.

Mapped initial function data are lifted onto prescribed and dependent traces
before starting native histories. The physical surface-L2 boundary projection
makes prescribed velocities have zero net mass flux while retaining side,
edge and corner traces. Its reported `boundary_flux_correction` is the RMS
top-surface velocity change; it does not add a volume mass source. Accepted
stages retain the original weak continuity rows and negligible gauge multiplier.
`RANSWallOwnershipShift` reassigns physical cells by
`(geometry ID+shift)%ranks` for partition/restart checks and is excluded from
the physics signature.

Mapped output uses `.chart.fld` with explicit coordinate/map metadata, together
with per-field `.vtu` files containing actual Cartesian points. Load the VTU
files for physical visualization. The JSON adds map parameters, roughness
force/work, the complete resolved-plus-k stage energy budget, full unmasked
final-stage loads, both regional stage loads and the associated constraints.
The energy budget uses the native combined stage history, which already
contains the explicit prescribed-body contribution.

The [M5 report](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-m5.md)
records geometry certificates, metric identities, normal traction and epsilon
constraint reactions, fixed-geometry refinement, partition/MPI and restart
evidence. The original model's stability and general inner-k positivity
limitations remain in force.


`roughness_budget` stores integrated Cartesian Fx,Fy,Fz, integrated U·F and
maximum |n·F|. `energy_budget` stores the unmasked velocity-plus-constant-k weak
test, stage storage, stress work minus k sources, convective work, prescribed
and roughness-force work, a reserved zero, their direct sum, weak-minus-direct
closure, and stored integral of |U|²/2+k. The energy audit includes essential
row reactions and does not assert unconditional physical decay.
