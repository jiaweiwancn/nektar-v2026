# Milestone 3 integration fixtures

These are small full-3D solver inputs for inlet integration, independent of the
Milestone 1 experimental data. They do not model the downstream experimental
decay. Numerical conclusions must use the Linux build and solver environment
described in the repository's `AGENTS.md`.

`HITInflowMeanOnly3D` and `HITInflowManufactured3D` are registered solver
regressions. `HITBoundaryUnitTests` exercises the actual boundary manager,
strict configuration and geometry guards, constrained projection, repeatability,
normal-velocity acceleration and projected-flux spatial convergence. The full
temporal-order study below is a separate numerical acceptance check.

`generate_m3_cases.py` uses only the Python standard library. It regenerates six
complete XML sessions, `m3-narrow-band.dat`, and a small independent oracle,
`m3-manufactured-reference.json`. The committed XML files can be loaded directly
by `SessionReader` and `ContField`; no mesh generator is needed. They contain four
affine hexes, two in x and two in y. Every input has the sparse inlet region
`REF="7"`, outlet region 11, and reciprocal periodic pairs 21/22 in y and 31/32
in z, for all velocity components and pressure. Pressure is the last variable.

| File | Purpose | Expansion modes | Timestep / steps |
| --- | --- | --- | --- |
| `m3-mean-only.xml` | Zero-RMS inlet and HOutflow; exact uniform `(1,0,0)`, `p=0` | 6 | 0.001 / 20 |
| `m3-stochastic.xml` | Nonzero von Karman inlet, mean initial field, HOutflow | 8 | 0.001 / 20 |
| `m3-manufactured-dt0.01000.xml` | Seed-matched manufactured flow | 12 | 0.01 / 10 |
| `m3-manufactured-dt0.00500.xml` | Same flow and final time | 12 | 0.005 / 20 |
| `m3-manufactured-dt0.00250.xml` | Same flow and final time | 12 | 0.0025 / 40 |
| `m3-manufactured-dt0.00125.xml` | Same flow and final time | 12 | 0.00125 / 80 |

The smoke cases are deliberately short. The stochastic case has periods
`Lx=Ly=Lz=2*pi`, a four-point synthesis lattice in each direction, and cutoff
1.55, retaining modes with nonzero kz as well as kx/ky. It has no `ExactSolution`;
the solver's printed error values are norms relative to zero. Their finite,
nonzero values establish a smoke run, not accuracy or experimental agreement.

## Manufactured solution

The manufactured cases have synthetic periods `Lx=Ly=2*pi`, `Lz=pi/2`, counts
`4,4,4`, target ensemble RMS 0.05, and seed 20261028. The linear table is zero
below 1.1, positive at `sqrt(2)`, and zero above 1.6. With cutoff 1.55 and taper
start 1.5, only the independent modes `(1,-1,0)` and `(1,1,0)` survive; the
negative conjugate partners are included by taking twice the real part. The
anisotropic two-mode support is chosen for a short analytical integration test.
It is insufficient for a physical HIT ensemble or spectrum assessment.

The generator independently implements the frozen `SplitMix64BoxMuller52V1`
counter, uniform conversion, complex Gaussian, and transverse cross-product
coefficient rule. Both modes have the same scalar spectral weight, so their
target-RMS amplitude factor simplifies to `sqrt(3*RMS^2/8)`. The reference JSON
records the integer counter words, coefficients, covariance and sample values.
It does not claim bitwise agreement between different libm implementations.

Let `v` denote the finite fluctuating Fourier field and `U=1`. The exact flow is

```
u(x,y,z,t) = U e_x + v(x-U*t,y,z),    p = 0.
```

Each mode is perpendicular to its wavevector, so this smooth field is
analytically divergence-free, and `d_t v + U d_x v = 0`. In the convention
`d_t u + u.grad(u) = -grad(p) + nu Laplacian(u) + f`, the body force is

```
f_i = sum_j v_j d_j v_i + nu sum_modes |k|^2 2 Re(a_i exp(i k.(x-U*t))).
```

The initial and exact solution functions contain those same seed-defined
coefficients. The HITInflow boundary supplies the exact incoming trace at x=0.
At x=1, the manufactured case uses the exact time-dependent normal derivative
of each velocity component and `p=0`. It therefore tests inlet integration and
the existing high-order pressure condition without introducing a modelled open
outlet error. The separate HOutflow cases exercise the intended open outlet.

Temporal acceptance uses the raw final-time velocity errors:

```
e(dt) = sqrt(L2(u)^2 + L2(v)^2 + L2(w)^2),
order = log(e(dt)/e(dt/2)) / log(2).
```

The default checker requires every adjacent order to be at least 1.8, and the
coarsest vector error to be at most `5e-5`. It reports individual velocity and
pressure L2/Linf errors as well. It does not subtract a fine-grid result or
relax a threshold when a test fails. The default mean-only threshold is `1e-10`
for all four variables' L2 and Linf errors. The stochastic check requires finite
bounded norms and nonzero transverse responses; it is intentionally weaker.

Twelve modes on the two y elements are intended to put spatial approximation
error below the temporal error. Confirm that on Linux before accepting an order
study: regenerate into a separate directory with `--mms-modes 14`, rerun at
least the finest timestep, and compare raw errors. Startup treatment of the
existing high-order pressure boundary can affect convergence. A failed order
must be investigated in the boundary/time-integration code and reported;
changing the exact inlet, suppressing the first errors, or reducing the required
order is not an acceptable substitute. This file records the test design, not
results from such a study.

## Linux run and analysis

Use a separate run directory outside the source checkout. Run sequentially
with one process and one thread initially; MPI comparisons can then reuse the
same inputs with an explicitly bounded rank count. Preserve complete logs and
field output in the run directory. For example, in a Linux shell:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
hit_fixture_dir='/media/jay/新加卷/nektar-v2026/solvers/IncNavierStokesSolver/Tests/HITInflow'
hit_run_dir='/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-3/cases'
python3 "$hit_fixture_dir/generate_m3_cases.py" --output-dir "$hit_run_dir"
cd "$hit_run_dir"
IncNavierStokesSolver m3-mean-only.xml > m3-mean-only.log 2>&1
python3 "$hit_fixture_dir/check_m3_cases.py" mean-only m3-mean-only.log
IncNavierStokesSolver m3-stochastic.xml > m3-stochastic.log 2>&1
python3 "$hit_fixture_dir/check_m3_cases.py" stochastic m3-stochastic.log
IncNavierStokesSolver m3-manufactured-dt0.01000.xml > m3-manufactured-dt0.01000.log 2>&1
IncNavierStokesSolver m3-manufactured-dt0.00500.xml > m3-manufactured-dt0.00500.log 2>&1
IncNavierStokesSolver m3-manufactured-dt0.00250.xml > m3-manufactured-dt0.00250.log 2>&1
IncNavierStokesSolver m3-manufactured-dt0.00125.xml > m3-manufactured-dt0.00125.log 2>&1
python3 "$hit_fixture_dir/check_m3_cases.py" --json temporal-result.json temporal \
    --run m3-manufactured-dt0.01000.xml m3-manufactured-dt0.01000.log \
    --run m3-manufactured-dt0.00500.xml m3-manufactured-dt0.00500.log \
    --run m3-manufactured-dt0.00250.xml m3-manufactured-dt0.00250.log \
    --run m3-manufactured-dt0.00125.xml m3-manufactured-dt0.00125.log
```

Record the tested source revision, build command, solver commands, MPI ranks,
thread settings, and any concurrent work against the 64-worker server cap.
For tabulated input, launch from the generated input directory so the relative
spectrum path resolves there. The scripts leave the original Milestone 1 data
and manifest unchanged.
