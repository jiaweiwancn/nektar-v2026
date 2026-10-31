# Milestone 5 focused HIT regressions

M5 delivered eleven focused HIT regressions when MPI and FFTW are enabled:
eight added in M5 and three mean/manufactured regressions from M3/M4. M6 added
three further checks, so the current configuration has fourteen solver entries.
The table below records the M5 subset; the current execution workflow is in
[the test/data guide](README.md). All use native executables and standard Tester
metrics. Running these solver tests does not require Python or optional unit-test builds.
The fixture generators use only the Python standard library.

| Test suffix after `IncNavierStokesSolver_` | Behavior checked |
| --- | --- |
| `HITInflowMeanOnly3D`, `HITInflowHomogeneous1DMean` | Exact zero-RMS mean limit |
| `HITInflowManufactured3D` | Two oblique modes, exact body force and time-dependent inlet |
| `HITInflow_m5_single_mode` | One independently specified transverse pair with known waveform |
| `HITInflowSeededH1D` | Seeded multimode von Karman velocity/pressure norms |
| `HITInflowTabulatedH1D` | Tabulated shell-spectrum velocity/pressure norms |
| `HITInflowProjectionH1D` | Raw/projected trace, nonzero packed z components and repeatability |
| `HITInflowProjectionH1D_FFTW` | The same native projection checks with FFTW |
| `HITInflowRestartH1D` | Preferred Checkpoint-filter restart versus uninterrupted flow |
| `HITInflowMPIH1D` | Serial versus eight-rank mixed mesh/Fourier flow, including inlet-less ranks |
| `HITInflowDivergenceFlux3D` | Manufactured bulk divergence, net flux and Gauss consistency |

FFTW and MPI cases are registered only when those features are enabled. The
MPI case declares eight CTest processors; its segments are sequential, and its
MPI command is first with an explicit `-n 8`. This ordering avoids the existing
Tester's serial-first/MPI-later launcher limitation. Its final serial native
comparator checks every field partition. Restart and MPI metrics require
matching identity/time, nonzero imported velocity components, and bounded
physical field differences. A command's zero exit status alone is insufficient.

The exact one-pair fixture uses a full-3D unit cube and synthesis periods
`(2*pi,1,1)`, counts `(4,4,4)`, seed 20261028, and cutoff 1.55. It retains only
canonical `(1,0,0)` (mode ID 58). The solution is
`(1, 2*Re(a_y*exp(i*(x-t))), 2*Re(a_z*exp(i*(x-t))))`, with zero pressure,
body force `(0,nu*v,nu*w)`, exact initial data and analytical outlet derivatives.
The golden coefficient/time/derivative values and fixed tolerances are in
`m5-single-mode-reference.json`. This deliberately unequal-period one-pair test
is manufactured, and does not claim isotropic ensemble statistics. It does not
use Nektar's unsupported `ModeType=SingleMode` expansion.

`TestHITSingleMode.cpp` adds four `HITCore` cases. They check normalized modal
transversality and `abs(trace(gradient))/FrobeniusNorm(gradient)` over nine
length/velocity rescalings, dimensionless finite-difference derivatives and
frozen convection, explicit zero normalization, and the independent golden
waveform. No dimensionful `max(1,scale)` floor is used. The core suite has 23
cases, alongside seven boundary and nine restart cases: 39 total.

`m5-acceptance.json` preserves the inherited/analytical limits. The separately
calibrated `m5-flow-acceptance.json` records the manufactured p/dt pilot and the
bulk-flow bounds. Bulk divergence is normalized as
`cbrt(volume)*L2(div u)/sqrt(sum_i L2(u_i)^2)`. Net-flux error is normalized by
`0.5*(abs(Qin)+abs(Qout))`; the Gauss residual compares the integral of bulk
divergence with the extracted-solution boundary flux. `HITFlowDiagnostics`
reads full native fields with coverage checks and integrates physical xy planes
with `LZ/N` weights for Homogeneous1D. Its flux uses solution traces, including
at the outlet, rather than the boundary-condition value arrays.

Build the solver, Tester and native tools with `NEKTAR_BUILD_TESTS=ON`.
`NEKTAR_BUILD_UNIT_TESTS=ON` additionally builds the three focused unit suites.
`HITBoundaryProbe`, `HITFieldCompare` and `HITFlowDiagnostics` now belong to the
solver test build and share the solver's binary directory; their source files
remain under `library/UnitTests/SolverUtils/`. Both field tools accept native
session options before their trailing field-file arguments.

For the authoritative Linux environment, activate the installed environment and
set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. Use `ctest -N` to inspect named
registrations. Run with explicit bounded concurrency; an eight-rank test must
be included in the combined 64-worker limit. The M5 validation exports the
registered Tester commands to separate working directories below the Linux run
root, preserving their arguments while keeping generated fields/logs outside
the checkout and home build tree. The closure report/receipt record that exact
runner and the selected legacy inflow, outflow, Homogeneous1D, transform and
moving-frame tests. These regressions establish numerical behavior on small
inputs. The completed limited-band production comparison and its limitations
are recorded in [Milestone 7](../../../../docs/plans/hit-inflow-milestone-7.md).
