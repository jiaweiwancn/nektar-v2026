# HIT inflow on a full 3D mesh

This complete session uses four hexahedra, degree 7, and all three velocity
components. A frozen von Karman realization is advected through the x-minimum
inlet at Ubar=1 m/s. The domain is 1 by 2*pi by 2*pi m. Both y and z faces are
periodic; the outlet uses `HOutflow` velocity-Neumann/pressure-Dirichlet conditions
and the inlet pressure uses `H`.

This is a small configuration demonstration, with nu=0.01 m2/s and a coarse
1.55 rad/m radial cutoff. It is not the grid-turbulence experiment or a
statistically converged turbulence calculation. The synthesis lattice is
(4,4,4), with a 2*pi m period in each direction, seed 20261028, and a cosine
energy taper starting at 1.5 rad/m. `TargetRMS=0.05 m/s` specifies expected
K=3*(0.05)^2/2; each finite realization can differ. SVV is off and spectral/hp
dealiasing is on. The initial volume contains the mean flow, so startup is part
of this 20-step example.

## Run in a fresh directory

Activate an installation containing `IncNavierStokesSolver`. On the project
Linux server:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
hit_src=/media/jay/新加卷/nektar-v2026
hit_run=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/examples/full3d-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$hit_run/tmp"
cp "$hit_src/solvers/IncNavierStokesSolver/Examples/HITInflow3D/HITInflow3D.xml" "$hit_run/"
cd "$hit_run"
export TMPDIR="$hit_run/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
IncNavierStokesSolver HITInflow3D.xml > solver.log 2>&1
```

This command uses one MPI rank. A two-rank run is also supported; use
`mpirun --bind-to none -n 2 IncNavierStokesSolver HITInflow3D.xml` in a separate
fresh directory. Account for all other work before allocating ranks: the
project Linux server permits at most 64 workers in total.

The solver should finish 20 steps at t=0.02 s. It writes `HITInflow3D.fld`,
checkpoint files/directories, an optimisation `.opt` file and `history.his`.
There are 12 probes and 21 equal-time frames, including the initial state.
The log reports input identity and expected/realized synthetic energy; generic
`L 2 error`/`L inf error` output is not an experimental accuracy measure.

With the optional test utilities built and installed:

```bash
HITFlowDiagnostics HITInflow3D.xml HITInflow3D.fld > diagnostics.log 2>&1
```

This saves `HITInflow3D.flow.json` with native divergence and integrated flux
checks. Those checks do not establish turbulence convergence. Plot or convert
the field using the normal `FieldConvert` workflow. For a restart, retain all
u/v/w/p fields and the embedded HIT identity; do not reset its absolute time.
The [restart regressions](../../Tests/HITInflow/README-m4.md) give complete cases.

## Configuration and experimental data

The user-guide subsection `sec:hit-inflow` documents the complete version-1
interface. The [HIT guide and data index](../../Tests/HITInflow/README.md)
explains shell versus longitudinal spectra, data conversion and tests.
The [Quasi-3D example](../HITInflowQuasi3D/README.md) uses a tabulated toy shell
spectrum with a Fourier direction. Both use the same synthesis conventions.
The [M7 report](../../../../docs/plans/hit-inflow-milestone-7.md) records the
limited-band experimental comparison and its failed precision criteria.

An initial legacy `_0.chk` at t=0 can precede HIT metadata setup and cannot be
used as a verified HIT restart when identity metadata is absent. The evolved
checkpoints at t=0.01/0.02 s contain the identity; retain all four fields and
check the stored time before restarting.
