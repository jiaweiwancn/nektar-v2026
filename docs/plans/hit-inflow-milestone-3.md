# HIT inflow Milestone 3 implementation and closure

Milestone 3 is closed on 2026-10-28. The full three-dimensional HIT inlet is
implemented, built and installed on Linux. All 26 unit cases and eight selected
solver regressions pass. The manufactured inlet retains second-order velocity
accuracy with the existing IMEX2 pressure startup policy.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | implementation and numerical verification |
| Origin date | 2026-10-28 |
| Version label | hit_inflow_full3d_v1 |
| Verification status | VERIFIED full-3D inlet, temporal/spatial checks, regressions and preliminary MPI setup on Linux |
| Final implementation revision | `961362b2d3c28758aac1ba68ca81a67b9f91bcaf` |
| Temporal study revision | `98fb4c8` |
| Spatial refinement and initial MPI checks | `974a973f6e8ab109ad32cef34bebe5b3785048f1` |
| Branch | `codex/hit-inflow-generator` |

The [machine-readable Linux receipt](hit-inflow-milestone-3-validation.json)
records full revisions, source and log hashes, commands, acceptance criteria and
results. Between the temporal-study revision and the final implementation, the
runtime changes only route agreed configuration/geometry failures through
Nektar's normal error handler. The generator, velocity projection, update
scheduling, pressure code and valid-input numerical paths are unchanged.

## Delivered implementation

Paths in this report are relative to the repository root.

| Component | Files | Behavior |
| --- | --- | --- |
| Vector inlet | `solvers/IncNavierStokesSolver/BoundaryConditions/HITInflow.h/.cpp` | Holds all three velocity expansions, caches quadrature coordinates, samples the immutable generator and applies constrained boundary transforms |
| Runtime configuration | `BoundaryConditions/HITInflowConfig.h/.cpp` in that solver | Strict version-1 XML, constant mean expressions, solver/integrator guards, root-only spectrum-file reading, broadcast of resolved settings and exact file bytes, SHA256 and diagnostics |
| Domain checks | `BoundaryConditions/HITInflowGeometry.h/.cpp` | Global REF to local-array mapping; collective geometry, normal, area, periodic-pair and outlet checks |
| Scheduling | `IncBaseCondition.h`, `IncBoundaryConditions.h/.cpp`, `EquationSystems/IncNavierStokes.cpp` | Separate prescribed-velocity callback at the end of `SetBoundaryConditions(time)` |
| Permanent checks | `library/UnitTests/SolverUtils/TestHITBoundary.cpp`, solver `Tests/HITInflowMeanOnly3D.tst` and `HITInflowManufactured3D.tst` | Actual boundary-object unit tests plus mean-only and manufactured solver regressions |
| Examples and numerical study | `solvers/IncNavierStokesSolver/Tests/HITInflow/` | Six complete sessions, a shell spectrum, an independent seeded oracle, fixture generator, raw-error checker and run instructions |
| User documentation | `docs/user-guide/solvers/incompressible-ns.tex`, HIT inlet subsection | Supported setup, XML, spectrum format, mean, record duration and projection limitations |

The boundary manager loads configuration on every session rank before its local
factory loop. It constructs the boundary-region description collectively and
maps ascending global REF IDs to the solver's local arrays, including periodic
regions. The tests deliberately use sparse inlet REF 7, whose local index is 0
in the serial fixture. A local inlet constructor starts no session collective.
One immutable realization is shared by the local inlet objects; evaluating a
time does not advance random state.

The hook supplies the initial trace before Dirichlet enforcement and the
next-step trace before the implicit velocity solve. Pressure extrapolation then
reads the projected boundary coefficients through its existing normal-velocity
history path. There is no extra `rhsTime+dt` inlet update. The ordinary boundary
manager update remains after pressure evaluation; the new inlet inherits its
empty implementation. Existing wall and moving-frame pressure updates retain
their scheduling.

The native pressure acceleration startup is also unchanged: its first two
pressure calls suppress the BDF acceleration term. The final-time manufactured
study below establishes second-order velocity convergence with that policy;
it does not assert exact pressure at the first startup step or restoration of
all multistep history on restart.

## Supported configuration and boundaries

The version-1 contract from [Milestone 1](hit-inflow-milestone-1.md) is now
implemented for full 3D, with these explicit staging limits:

- Use `VelocityCorrectionScheme`, `UnsteadyNavierStokes`, nonlinear advection,
  continuous `u,v,w,p`, `Standard` extrapolation, fixed-step IMEX order 2 with an
  empty variant, and no flow-rate/adaptive-step controller.
- Use one complete, axis-aligned boundary region on each side of a rectangular
  box. The inlet lies at `x=ORIGIN.X`, on the box's x minimum, with outward normal
  `(-1,0,0)`. All three velocity conditions are `D/HITInflow`; pressure is
  `N/H` with XML value zero.
- Pair both y faces and both z faces periodically for all four variables. Pairs
  must be reciprocal and translational, and match `LY` and `LZ`. The synthetic
  streamwise period need not equal the open domain length.
- At the x-maximum outlet, use velocity Neumann conditions, or the existing
  `HOutflow` Robin form, with pressure Dirichlet or `HOutflow` Robin conditions.
  The smoke fixtures use the standard Neumann/Dirichlet `HOutflow` combination.
- Ordinary `Body` momentum forcing is permitted for manufactured solutions.
  The convection speed remains the constant inlet u expression. Other forcing
  types, moving frames, mappings, extra variables and unsupported schemes fail
  before mode allocation.

`Homogeneous1D` is deliberately rejected until Milestone 4 implements and tests
the physical/Fourier boundary transform path. The preliminary full-3D MPI
checks here do not close the broader Milestone 4 decomposition and restart
requirements.

The parser checks complete job duration against the declared frozen record,
cutoff/taper and synthesis Nyquist limits, and the temporal sampling limit.
The seed is parsed as unsigned decimal uint64, separately from floating-point
session parameters. Unknown/duplicate custom elements and attributes,
nonconstant or random scalar expressions and inconsistent mean specifications
fail. Table input uses the solver working directory and the strict SI
`Shell3D` format. Startup reports resolved settings, mode memory before
allocation, expected/realized energy and the exact-byte digest for table input.

## Linux verification

All compilation, solver runs and numerical acceptance checks used
`jay@jay-Super-Server`, source `/media/jay/新加卷/nektar-v2026`, Release MPI build
`/home/jay/build/nektar-v2026-mpi`, and installation
`/home/jay/.local/opt/nektar-v2026-mpi`. Compiler: GCC 11.4, C++17, `-O3 -DNDEBUG`.
MPI is enabled, FFTW is disabled, and separate solver libraries are disabled.

Compilation used 16 workers. Unit tests and regression runners used concurrency
1. Every numerical process used `OMP_NUM_THREADS=1` and
`OPENBLAS_NUM_THREADS=1`; MPI smoke jobs used 4 ranks. The largest concurrent
solver allocation was six workers, and no build overlapped a solver run.
All solver data, logs and heavy temporary files are under
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-3/`. Test inputs and
collection-autotuning output also use temporary storage outside the source
checkout. The two checkouts' existing files were preserved.

### Core and boundary unit checks

The 19 existing core cases and 7 new boundary cases pass. Boundary coverage
includes three SHA256 known answers (empty, `abc`, one million `a` bytes), 13
invalid-configuration subcases, two inconsistent-geometry subcases, optional
configuration in a non-HIT session, exact mean-only projection, all three
nonzero velocity components, repeated/backward-in-time sampling, untouched
pressure coefficients, and nonzero prescribed normal acceleration checked
against analytical derivatives.

For the cubic nine-mode realization, the integrated projected streamwise
fluctuation over the inlet changes from `-7.20296e-8` to `7.82715e-10 m^3/s`
when expansion modes increase from 12 to 14. The refined result is below the
`2e-8 m^3/s` acceptance bound and below one tenth of the coarse magnitude.
The realization's mean is never adjusted to enforce a quadrature flux.

This distinguishes the analytically solenoidal Fourier field from its finite
polynomial trace. Bulk discrete divergence and resolved statistical convergence
remain later validation tasks.

### Solver accuracy and outlet checks

The mean-only HOutflow case finishes 20 steps to `t=0.02 s`. The largest L2/Linf
error over `u,v,w,p` is `5.30622e-11`, below `1e-10`. The stochastic von Karman
case also finishes, with finite velocity/pressure norms and nonzero transverse
velocities. Its mean-only initial field has startup transients; this short run
establishes outlet compatibility, not stationary HIT statistics.

The manufactured test convects two oblique, seed-defined transverse modes with
mean speed `1 m/s`. It uses exact initial values, analytical body forcing,
the generated inlet, and an exact time-dependent velocity-Neumann/zero-pressure
outlet. The separate HOutflow cases test the practical open boundary.

All manufactured cases finish at `t=0.1 s`. The table uses raw solver-reported
L2 errors, combined as `sqrt(L2(u)^2+L2(v)^2+L2(w)^2)`, without subtracting a
reference run:

| dt (s) | Velocity vector L2 | Order from previous row | Pressure L2 |
| --- | --- | --- | --- |
| 0.01 | 2.61807774e-5 | — | 6.49521e-6 |
| 0.005 | 6.14171048e-6 | 2.0918 | 7.33071e-7 |
| 0.0025 | 1.53688870e-6 | 1.9986 | 1.19742e-7 |
| 0.00125 | 3.85425010e-7 | 1.9955 | 2.49858e-8 |

Every adjacent velocity order exceeds 1.8 and the coarse error is below
`5e-5`. Increasing expansion modes from 12 to 14 at the finest timestep changes
the velocity vector error to `3.86002633e-7`, a relative change of 0.1499%, below
the 1% acceptance criterion. Pressure L2 is `1.63706e-8` in that refined run.
This supports the temporal study without a dominant spatial-error floor.

### Regressions and collective setup

The two permanent HIT solver regressions pass, together with these six existing
ones, using their unchanged metrics:

- `Womersley_PipeFlow`
- `Tet_channel_m3_ConOBC`
- `Tet_Kovasnay_HOBC`
- `MovingRefFrame_translation`
- `MovingRefFrame_Rot_SimpleDomain`
- `ChanFlow_3DH1D_bcsfromfiles`

Four-rank full-3D smoke runs pass with both von Karman and tabulated spectra.
Partitions 0 and 3 have no inlet faces; rank 0 still loads and distributes the
configuration and table. The table digest matches an independent SHA256 of
its bytes. A missing table is rejected collectively before stepping with exit
status 1, without an uncaught exception or timeout. These are setup and finite-response checks; exact
cross-decomposition field reproducibility remains in Milestone 4.

## Reproducing the checks

Activate the Linux installation and set both thread counts to 1. Build and
install with:

```bash
cmake --build /home/jay/build/nektar-v2026-mpi \
    --target IncNavierStokesSolver SolverUtilsUnitTests HITBoundaryUnitTests \
    --parallel 16
ctest --test-dir /home/jay/build/nektar-v2026-mpi \
    -R '^(HITCoreUnitTests|HITBoundaryUnitTests)$' --parallel 1 --output-on-failure
cmake --install /home/jay/build/nektar-v2026-mpi
```

Use a separate run directory and follow
[`Tests/HITInflow/README-m3.md`](../../solvers/IncNavierStokesSolver/Tests/HITInflow/README-m3.md)
for the six serial inputs and raw-error checker. Regenerate into a different
directory with `--mms-modes 14` for the spatial check. The receipt records exact
commands, including `Tester --executable ...` runs from external storage and
bounded four-rank MPI commands with stdin redirected from `/dev/null`.

Milestone 4 remains open for Homogeneous1D sampling/transforms, broader MPI
reproducibility and restart work. No experimental downstream-decay comparison
or production-resolution HIT validation is claimed by this milestone.
