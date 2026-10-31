# Milestone 4 Homogeneous1D and MPI fixtures

These are small integration inputs. They use the same immutable seeded HIT
realization in full 3D and Homogeneous1D; they do not model the experimental
downstream-decay case. All compilation, solver runs and numerical conclusions
must use the authoritative Linux environment in `AGENTS.md`.

`generate_m4_cases.py` uses the Python standard library and the existing M3
fixture generator. It writes 13 complete sessions and a shell-spectrum table.
It leaves every M1 experimental input and M3 fixture unchanged. Default inputs
use four xy quadrilaterals or their four-hex extrusion, expansion modes 8,
eight homogeneous Fourier planes, `dt=0.001`, and 20 uninterrupted steps to
`t=0.02`. The box spans x in `[0,1]`, y and z in `[0,2*pi]`. Both representations
use synthesis periods `(2*pi,2*pi,2*pi)`, counts `(4,4,4)`, seed 20261028,
cutoff 1.55, taper start 1.5, and ensemble target RMS 0.05.

| Input | Purpose |
| --- | --- |
| `m4-h1d-von-karman.xml` | Homogeneous1D baseline, matrix Fourier transforms |
| `m4-full3d-von-karman.xml` | Same physical realization on a full 3D mesh |
| `m4-h1d-mean-only.xml`, `m4-full3d-mean-only.xml` | Zero-RMS exact uniform-flow limit |
| `m4-h1d-tabulated.xml`, `m4-full3d-tabulated.xml` | Tabulated shell-spectrum path |
| `m4-h1d-fftw.xml` | Same H1D von Karman case with `USEFFT=FFTW` |
| `m4-h1d-segment.xml`, `m4-full3d-segment.xml` | First half of the record, checkpoint at `t=0.01` |
| `m4-h1d-restart.xml`, `m4-full3d-restart.xml` | Read the corresponding segment checkpoint and advance the remaining 10 steps |
| `m4-h1d-checkpoint-filter.xml`, `m4-h1d-filter-restart.xml` | Exercise metadata and restored time through the preferred Checkpoint filter |

Stochastic cases use `IterativeStaticCond` and
`IterativeSolverTolerance=1e-14` to make serial/MPI comparisons meaningful at a
common requested solver tolerance. The exact uniform-flow cases use
`DirectStaticCond` to isolate DC preservation from iterative pressure-solve
noise; their error bound remains `1e-10` for every component and pressure.
The predeclared numerical acceptance bounds are in `m4-acceptance.json`.
Sparse inlet REF 7 prescribes all three velocity components and high-order
pressure data. Outlet REF 11 uses the existing Neumann/Dirichlet HOutflow
condition. The y pair is 21/22; full 3D also has explicit z pair 31/32, whereas
Homogeneous1D uses the Fourier period `LZ`.

The spectrum file `m4-shell-spectrum.dat` is an exact copy of the existing
manufactured integration-test table `m3-narrow-band.dat`. It is not experimental
data. On the cubic synthesis lattice its nonzero band retains six independent
oblique modes, including nonzero z wavenumbers. The von Karman case retains
nine independent modes. A finite small lattice is useful for transform and MPI
checks, not for a resolved physical HIT spectrum.

Generate a separate Linux run directory outside the checkout:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
hit_m4_fixture='/media/jay/新加卷/nektar-v2026/solvers/IncNavierStokesSolver/Tests/HITInflow'
hit_m4_run='/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-4/cases'
python3 "$hit_m4_fixture/generate_m4_cases.py" --output-dir "$hit_m4_run"
cd "$hit_m4_run"
IncNavierStokesSolver m4-h1d-von-karman.xml > h1d-serial.log 2>&1
mpirun -np 4 IncNavierStokesSolver --npz 1 m4-h1d-von-karman.xml < /dev/null > h1d-mesh4.log 2>&1
mpirun -np 4 IncNavierStokesSolver --npz 4 m4-h1d-von-karman.xml < /dev/null > h1d-homogeneous4.log 2>&1
mpirun -np 4 IncNavierStokesSolver --npz 2 m4-h1d-von-karman.xml < /dev/null > h1d-mixed4.log 2>&1
```

Run these commands in separate output directories if preserving field files or
partition XML is necessary: repeated runs with the same session stem otherwise
reuse output names. Check the actual worker allocation before concurrent runs;
the repository's combined Linux limit is 64 workers.

For a restart comparison, run the uninterrupted baseline and the segment, then
the restart from the same directory. The segment uses legacy
`IO_CheckSteps=10` and creates `<segment-stem>_0.chk` and
`<segment-stem>_1.chk`; the latter is the `t=0.01` restart input. Both segment
and restart retain the full frozen-record duration 0.02 and time origin zero.
The restored field supplies the actual starting time. The generator boundary
must be reconstructed at that time, while the solver's native IMEX history is
reinitialized; an evolved restart is therefore a tolerance/convergence check,
not a promise of bitwise uninterrupted flow.

Use `--modes 10` for polynomial refinement and `--homogeneous-modes 16` for a
Fourier refinement input. Use `--dt 0.0005 --steps 40` to halve the timestep at
the same final time; the generated segment/restart pair splits that record in
half. FFTW inputs require a build with `NEKTAR_USE_FFTW=ON`. This file describes
test inputs and reproducible commands; numerical pass/fail evidence belongs in
the Milestone 4 closure report and validation receipt.

The native `HITBoundaryProbe` accepts the usual session and MPI partition
arguments. Run it from the generated-input directory; it writes one CSV per
rank and a `.probe.json` summary beside the copied session. It samples all
components at five times, including repeated/backward evaluations, and inverts
the actual constrained projection to physical space. Use `check_m4_cases.py
probe` and `samples` to verify packed nonzero Fourier content and compare
common coordinates across representations, MPI partitions and transform
backends.

`HITFieldCompare session.xml first.fld second.fld` runs with one process and
imports every partition of both field directories. It checks representation,
period, polynomial modes and complete element/Fourier-plane coverage before
reporting physical L2/Linf differences at full printed precision. Use the
checker's `fields` command for same-representation MPI/backend comparisons;
restart comparisons use the separate bounds and timestep-refinement criterion
in the acceptance file. Raw solver norms can be checked with `solver` and
compared between full 3D and Homogeneous1D with `statistics`.
