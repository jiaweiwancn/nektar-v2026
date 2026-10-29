# HIT inflow Milestone 5 closure: focused regression coverage

Milestone 5 is **CLOSED on 2026-10-29**. Eight durable solver regressions and
four core cases were added. All 39 focused unit cases and 19 selected solver
regressions pass on Linux, including the full solver regression set with
optional unit tests disabled. Bulk-divergence and flux limits are calibrated
from manufactured-solution refinement and checked by a fresh solve.

This closes Milestone 5 of the [implementation plan](hit-inflow-implementation-plan.md).
Production generator, boundary sampling and solver equations are unchanged
from Milestone 4. The work adds test inputs, checking tools, registration and
reproducibility evidence.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | implementation and numerical verification |
| Origin date | 2026-10-29 |
| Version label | hit_inflow_regressions_v1 |
| Verification status | VERIFIED for the stated Linux cases |
| Final implementation revision | `ce05c35edbf737eecc8e2aab1561ac57294199b1` |
| Unit/build candidate | `46e4da4` |
| Branch | `codex/hit-inflow-generator` |
| Numerical receipt | [Revisions, commands, source/log hashes and measurements](hit-inflow-milestone-5-validation.json) |

The final complete 19-test batch uses `ce05c35`. The 39 unit cases use
`46e4da4`; subsequent changes affect native test-helper argument/validation
logic and calibrated regression definitions, with no change to the core or
boundary/restart unit sources. The receipt records full revisions separately.

## Durable coverage

All cases use `ADD_NEKTAR_TEST`; MPI/FFTW-specific tests are conditional on
those features. The three existing HIT mean/manufactured tests remain.
The following eight are new:

| Test suffix after `IncNavierStokesSolver_` | Check |
| --- | --- |
| `HITInflow_m5_single_mode` | One transverse pair, independent waveform and exact manufactured solution |
| `HITInflowSeededH1D` | Seeded multimode von Karman field, including nonzero transverse components |
| `HITInflowTabulatedH1D` | Supplied three-dimensional shell spectrum |
| `HITInflowProjectionH1D` | Raw and constrained trace, packed nonzero z content, repeatability |
| `HITInflowProjectionH1D_FFTW` | The same projection checks with FFTW |
| `HITInflowRestartH1D` | Preferred Checkpoint-filter restart versus uninterrupted flow |
| `HITInflowMPIH1D` | Serial versus eight-rank mixed mesh/Fourier flow |
| `HITInflowDivergenceFlux3D` | Manufactured bulk divergence, net flux and Gauss consistency |

The [fixture instructions](../../solvers/IncNavierStokesSolver/Tests/HITInflow/README-m5.md)
list the complete eleven HIT regressions, input files, equations, bounds and
native-tool behavior. Most new cases reuse the small committed M3/M4 inputs;
M5 adds an exact one-pair session and an MPI session with a separate output stem.
The tabulated spectrum is manufactured test data, not the experimental spectrum.

### One independently specified transverse pair

The full-3D unit cube uses synthesis periods `(2*pi,1,1)`, counts `(4,4,4)`,
seed 20261028 and cutoff 1.55. Only canonical index `(1,0,0)`, mode ID 58,
survives. With the independently computed complex coefficients, the exact
solution is

`u=1`, `v=2*Re(a_y*exp(i*(x-t)))`, `w=2*Re(a_z*exp(i*(x-t)))`, `p=0`.

Body force `(0,nu*v,nu*w)`, exact initial conditions and analytical outlet
Neumann data maintain that frozen advected wave. The fixture uses 20 IMEX2
steps of `dt=0.005` to `t=0.1`. Its fixed L2/Linf limits are `1e-9` for u/p,
`2e-6` for v/w L2, and `3e-6` for v/w Linf. A Regex metric requires exactly
one independent mode. Coefficients, counter words and phase/derivative samples
are saved in `m5-single-mode-reference.json` and checked by a golden unit case.

The unequal-period lattice deliberately isolates one pair. It is not a physical
isotropic ensemble and does not use Nektar's unsupported `ModeType=SingleMode`.

### Core invariants without a dimensional floor

`TestHITSingleMode.cpp` adds four cases to the existing `HITCore` suite:

- Modal `abs(k dot a)/(|k|*sqrt(sum(abs(a_i)^2)))` and sampled
  `abs(trace(gradient))/FrobeniusNorm(gradient)` below `1e-12`.
- Dimensionless spatial/time finite-difference derivative checks below `1e-8`,
  analytic frozen-convection residual below `1e-12`, and convected-sample agreement.
- Explicit zero-field normalization: zero residual over zero scale is zero;
  a nonzero residual cannot be hidden by a zero denominator.
- Golden one-pair coefficients and analytic waveform/derivatives near
  floating-point precision, normalized by the specified RMS.

The first two cover all nine combinations of length scales `1e-6,1,1e6` and
velocity scales `1e-8,1,1e8`. They use no dimensionful `max(1,scale)` floor.
Together with existing coverage, 23 core, seven boundary and nine restart
cases pass: 39 total. Existing projected-inlet flux refinement from p12 to p14
and its refined absolute flux bound remain in the boundary suite.

## Bulk divergence and flux calibration

`HITFlowDiagnostics` imports every field partition and differentiates the
physical reconstructed velocity with native operators. It reports

- `D = L2(div u)`;
- `eta_D = cbrt(volume)*D/sqrt(sum_i L2(u_i)^2)`;
- `eta_Q = abs(Qin+Qout)/(0.5*(abs(Qin)+abs(Qout)))`;
- a Gauss residual comparing `integral(div u)` with `Qin+Qout`, using the same
  through-flow normalization.

Flux is integrated from solution traces extracted from the imported field.
Outlet boundary-condition arrays are not substituted for outlet velocity.
Homogeneous1D diagnostics integrate physical xy planes with Fourier weight
`LZ/N`; implicit periodic z boundaries contribute zero net flux. Complete-field
import and serial diagnostic execution avoid boundary-subset collectives.

The M5 tool post-processes the already validated M3 manufactured fields at
`t=0.1`; their original revisions and exact field/XML hashes are in the receipt.
No new solver accuracy claim is inferred from their filenames alone.

| Retained manufactured solution | D | eta_D | eta_Q |
| --- | --- | --- | --- |
| p12-dt010 | 4.7667410e-06 | 3.2475906e-06 | 3.5675049e-10 |
| p12-dt005 | 3.0316093e-07 | 2.0654421e-07 | 1.6813972e-12 |
| p12-dt0025 | 5.0946799e-08 | 3.4710167e-08 | 1.5759274e-12 |
| p12-dt00125 | 1.0140601e-08 | 6.9088142e-09 | 3.2612843e-13 |
| p14-dt00125 | 6.4026169e-09 | 4.3621172e-09 | 4.4977646e-13 |

Timestep refinement decreases bulk divergence; p14 also reduces the finest-step
divergence. Flux reaches a small floating-point/solver-error range rather than
remaining monotonically decreasing at its final values. The pilot's Gauss
residuals are below `5e-16`.

The following limits were frozen in `m5-flow-acceptance.json` before the fresh
regression solve: `D<=1e-5`, `eta_D<=1e-5`, `eta_Q<=1e-8`, and Gauss residual
`<=1e-11`. They apply to the stated p12, `dt=0.01` manufactured case and its
native solver settings; they are not universal production-HIT tolerances.

The fresh regression passes with `D=4.76674090e-6`, `eta_D=3.24759058e-6`,
`eta_Q=3.56733208e-10`, and Gauss residual `9.13e-17`. The inherited/analytical
criteria in `m5-acceptance.json` remain byte-identical to the initial M5 commit;
the calibrated flow criteria remain byte-identical to `ce05c35`.

### Short stochastic Homogeneous1D diagnostic

Both matrix and FFTW paths were also post-processed using the retained short,
mean-initialized M4 field. They give `eta_D` about `1.29449e-3`, differing by
approximately `2.23e-14`, net-flux error about `1.01e-7`, and Gauss residuals
below `1.5e-16`. This checks the H1D derivative/trace implementation while
recording a real numerical limitation of that transient fixture. It does not
meet the manufactured case's separate divergence limit and is not treated as
a converged decay dataset. Milestone 6 must assess startup, resolution and
domain effects before physical statistics are interpreted.

## Restart, MPI and checker integrity

The preferred Checkpoint-filter comparison preserves the immutable identity
and final time. Restarted native IMEX history is reinitialized; the test checks
bounded evolved error rather than bitwise equality.

| Comparison | Maximum velocity Linf difference | Pressure Linf difference | Bounds |
| --- | --- | --- | --- |
| Checkpoint-filter restart | `1.15300950e-5` | `1.48308996e-5` | `2e-4`, `2e-3` |
| Serial versus mixed eight-rank flow | `3.35547397e-12` | `5.47424232e-12` | `2e-9`, `2e-8` |

The MPI layout is four mesh partitions times two Fourier partitions. With four
xy elements and two inlet edges, it exercises ranks without inlet faces.
Both comparisons require identity/time agreement and three nonzero imported
velocity components, in addition to full field coverage and numeric bounds.

Three negative checker validations pass:

1. Equal malformed timestamps such as `0.02suffix` cause the actual Tester
   identity/time metric to fail.
2. Two otherwise matching, valid-header fields with zero payloads cause the
   actual nonzero-component metric to fail. Zero differences alone cannot pass.
3. A valid single partition supplied instead of the complete field directory
   is rejected by the importer's element/Fourier-plane coverage check.

These validations mutate only copies in experiment storage. Native field tools
continue to print measurements; Tester applies the committed acceptance metrics.

## Test infrastructure and integration correction

The native checking targets now belong to the solver's `NEKTAR_BUILD_TESTS`
build and share its binary directory, independently of `NEKTAR_BUILD_UNIT_TESTS`.
Their source files remain under `library/UnitTests/SolverUtils/`.
`HITFieldIO.h` factors the prior coverage-checked native importer into one
shared helper, retaining native UNIORDER and period-header compatibility.

The first restart/MPI trial reached the final checker but failed because Tester
inserts `--use-opt-file test.opt`, while the old helper accepted exactly four
arguments. Both field tools now accept native session options before their
trailing field-file arguments. The fixed trials and final complete batch pass;
no acceptance limit was relaxed to address the interface failure.

The MPI .tst starts with the explicitly eight-rank segment. This avoids the
existing Tester's unsafe serial-first/MPI-later prefix behavior; M5 does not
change the general runner. The final checker is last because sequential
segments overwrite prior output files. The MPI case declares eight CTest
processors. `NEKTAR_TEST_FORCEMPIEXEC=OFF` is recorded for the tested environment.

## Authoritative Linux results and reproduction

Source: `/media/jay/新加卷/nektar-v2026`; build:
`/home/jay/build/nektar-v2026-mpi`; installation:
`/home/jay/.local/opt/nektar-v2026-mpi`. Evidence is retained under
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-5/46e4da4`.
The directory label identifies the initial candidate, not every tested revision.
The receipt hashes 85 retained artifacts; source M3/M4 field hashes
used for post-processing are recorded separately.

All nineteen registered solver tests pass in a single CTest batch:

- Eleven HIT cases listed in the fixture instructions.
- `ChanFlow_3DH1D_bcsfromfiles` and `ChanFlow_3DH1D_FFT`.
- `KovaFlow_3DH1D_MVM_FFTW_Consistency`.
- `MovingRefFrame_translation` and `MovingRefFrame_Rot_SimpleDomain`.
- `ChanFlow_3DH1D_Parallel_mode1` and `ChanFlow_3DH1D_Parallel_mode2`.
- `KovaFlow_m8_short_HOBC`, the small existing high-order outflow regression.

For this batch, unit tests were disabled, all native regression tools built,
and their exact registered Tester argument vectors exported from CTest JSON to
external per-test working directories. Only working directory and explicit
scheduler properties were adapted for the repository storage/worker policy.
Nineteen tests passed with eight scheduler slots; the unit-test option was then
restored to ON. No optional Python package is needed by the registered tests.
The external orchestration uses only Python's standard library.

Representative recorded commands are:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cmake --build /home/jay/build/nektar-v2026-mpi \
  --target IncNavierStokesSolver SolverUtilsUnitTests HITBoundaryUnitTests \
           HITBoundaryProbe HITFieldCompare HITFlowDiagnostics Tester \
  --parallel 32
ctest --test-dir /home/jay/build/nektar-v2026-mpi \
  -R '^(HITCoreUnitTests|HITBoundaryUnitTests|HITRestartUnitTests)$' \
  --parallel 1 --output-on-failure
```

The receipt also records the unit-disabled configure/build, external CTest
command, restoration, helper paths/hashes, actual test arguments and logs.
The Release environment has MPI/FFTW enabled and HDF5 disabled. Compiler jobs
are capped at 32; numerical runs use one OpenMP/BLAS thread. Final CTest uses
eight worker slots, and concurrent diagnostic/regression work never exceeds
nine validation workers. Builds do not overlap solver execution, remaining
below the hard 72-worker limit.

## Remaining scope

Milestones 6-8 remain open. The tests freeze reproducible numerical behavior
and useful diagnostics on small cases. Production resolution, domain/record
length, startup flushing, uncertainty, downstream turbulence decay and the
Comte-Bellot--Corrsin comparison remain later work. HDF5 checkpoint behavior
is not exercised by this XML-based build.
