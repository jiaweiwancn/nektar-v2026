# Revised CLS/EP-AWF native flat wall stage

`SolverType=RANSEquilibriumWallSystem` selects the revised-note wall element.
The implementation uses native continuous spectral elements for velocity and
k across the full domain, outer-only epsilon, and separate inner/outer reduced
pressure spaces. It applies one global volume-mean pressure gauge. Mesh
coordinates `(x,y)` represent physical `(x,z)` and variables `(u,v)` represent
`(U,W)`. The outer CLS stress retains all three modeled normal components.

This M3 backend solves one supplied nonlinear stage in serial on affine,
axis-aligned, canonically ordered quadrilaterals. The mesh must have a complete
element boundary at the constant matching height H. Polynomial bases are
Modified_A, with equal degrees in both directions and degree at least two;
the supplied fixtures use velocity/scalar degree p and pressure degree p−2.
The dense coupled backend accepts at most 3,000 unknowns. Unsupported meshes
and MPI partitions are rejected. Native IMEX stepping, restart, three-dimensional
fields and curved geometry belong to subsequent milestones.

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
