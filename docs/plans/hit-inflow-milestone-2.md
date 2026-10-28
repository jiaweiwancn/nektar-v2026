# HIT inflow Milestone 2 implementation and closure

Milestone 2 is closed on 2026-10-28. Nektar++ now has reusable spectrum providers
and a deterministic three-component Fourier synthesis core implementing the
isotropic Batchelor tensor. Its independent Linux checks pass for normalization,
covariance, real-valued synthesis, divergence, analytical derivatives, and
Gaussian ensemble statistics.

The inlet boundary class, XML runtime parser, pressure update hook, distributed
sampling, and downstream-decay simulation remain in Milestones 3-7. The
Milestone 1 XML examples do not yet run an inlet in the solver.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | implementation and numerical verification |
| Origin date | 2026-10-28 |
| Version label | hit_inflow_core_v1 |
| Verification status | VERIFIED core build and focused numerical tests on Linux |
| Tested revision | `a6790a4a505a78e363759383478599efe2bd837b` |
| Branch | `codex/hit-inflow-generator` |

## Delivered code

Repository-relative paths:

| Component | Files | Behavior |
| --- | --- | --- |
| Spectrum providers | `library/SolverUtils/Turbulence/HITSpectrum.h/.cpp` | von Karman and tabulated SI shell spectra, validation, interpolation, explicit low/high extensions |
| Fourier synthesis | `library/SolverUtils/Turbulence/HITGenerator.h/.cpp` | Canonical half-lattice modes, versioned Gaussian coefficients, transverse factorization, finite-lattice normalization, physical sampling and derivatives |
| Core tests | `library/UnitTests/SolverUtils/` | New `SolverUtilsUnitTests` executable and focused CTest `HITCoreUnitTests` |
| Build registration | `library/SolverUtils/CMakeLists.txt`, `library/UnitTests/CMakeLists.txt` | Builds the core in SolverUtils, installs its public headers, and adds the unit-test subdirectory |

The API uses `NekDouble`, standard C++17 arrays, complex values and vectors. It
has no XML, MPI collective, boundary-expansion, or FFT dependency. The reference
evaluator sums the retained Fourier series directly; its cost grows with the
number of modes and sample points. Later sampling optimizations must preserve
this reference result.

## Native API and conventions

The public classes are in `Nektar::SolverUtils`. A small native example is:

```cpp
#include <SolverUtils/Turbulence/HITGenerator.h>

using namespace Nektar::SolverUtils;

auto spectrum = HITSpectrum::VonKarman(0.3, 1.7);
HITGeneratorConfig config;
config.seed = HITGenerator::ParseSeed("20261028");
config.normalization = HITNormalization::TargetRMS;
config.targetRMS = 0.222;

// Preflight scans the spectrum without allocating modes or drawing randomness.
const auto modeBytes = HITGenerator::EstimateModeBytes(config, spectrum);
HITGenerator generator(config, spectrum);
const auto sample = generator.Evaluate({0.2, 0.3, 0.4}, 0.1);
const auto statistics = generator.Diagnostics(1.5e-5);
```

This example uses the API's small default lattice and record duration. It is a
sampling example, not a resolved experimental configuration. The native API
normalizes the supplied spectrum when `TargetRMS` is selected; the future XML
parser separately enforces the Milestone 1 rule that this option must not also
specify an amplitude.

`HITSpectrum::FromFile` and `ReadTable` require the versioned `Shell3D` header
and SI angular-wavenumber convention. They reject a `Longitudinal1D` spectrum,
wrong units, malformed decimal tokens, unordered/duplicate wavenumbers, extra
columns, nonfinite values, and incompatible interpolation/tail settings.
Endpoints retain their exact input values. Log-domain formulas handle very
large or small wavenumbers without unnecessarily overflowing intermediate
powers or Gaussian-tail coefficients. Mathematical underflow can return zero;
unrepresentable finite-range results fail explicitly.

`HITGeneratorConfig` contains box periods/counts, phase origin, seed, positive
mean speed, phase-origin time, record duration, fixed sample timestep, cutoff,
energy taper, repeat policy, and normalization rule. Its optional
`maxModeBytes` is a native resource guard, independent of random realization
identity. It rejects an excessive mode allocation before reserving that vector.

`GetModes()` returns only canonical members of conjugate pairs. Sampling uses
twice the real part, giving the complete real field. Nyquist planes and the
zero vector are excluded. `Evaluate(point,time)` returns:

- `velocity`, including the constant mean in the x component;
- `gradient[i][j]`, the analytical derivative of velocity component i with
  respect to coordinate j;
- `timeDerivative`, the analytical derivative of the frozen convected field.

The derivative identity is `timeDerivative[i] = -Ubar*gradient[i][0]`. Evaluation
uses absolute time and does not consume RNG state. Reordering calls or
reconstructing a generator with the same settings gives the same sample on the
same numerical platform. Bitwise portability of libm transcendentals across
different platforms is not claimed.

`Diagnostics(viscosity)` reports fluctuating-field quantities, excluding the
mean kinetic energy: expected and realized kinetic energy, expected and realized
covariance tensors, and expected and realized viscous dissipation. The realized
values are exact Fourier/periodic-volume sums for that retained realization,
not finite probe-record estimates. Both use all conjugate partners consistently.
The reported `energyScale` is the actual scalar used in mode construction.

`TargetRMS` fixes the finite-lattice ensemble expectation
`expected(K)=3*RMS^2/2`; it does not force each realization's energy or component
variance to equal that value. A finite anisotropic lattice can retain component
differences. Exact zero target RMS or preserved zero amplitude uses the
mean-only path without mode allocation or Gaussian draws.

## Validation scope and independent checks

The 19 Boost cases cover:

1. Exact integer mixing, open-interval uniform endpoints, mode IDs and unsigned
   seed/counter wrap conventions.
2. Strict seed syntax and integer overflow/index validation.
3. von Karman values against the direct analytical formula and range handling.
4. Linear/log interpolation, exact endpoints, K4/zero low extensions and
   Gaussian/zero high extensions.
5. Invalid table shapes, values and enum policies.
6. Spectrum header/kind/units and decimal-token validation.
7. Original experimental table values, missing-point counts and rejection of a
   one-dimensional longitudinal spectrum as an inlet shell spectrum.
8. Canonical half-lattice content, mode transversality, and diagnostics against
   a separately summed full-lattice Batchelor tensor.
9. Ensemble normalization and seed-dependent realized energy.
10. Allocation-free mean-only fields and zero derivatives.
11. Real synthesis against an explicit complex sum of both conjugate members.
12. Analytical gradients/time derivatives against central differences, zero
    divergence, and invariance when following the mean convection.
13. Deterministic reconstruction, repeated and reordered calls, spatial and
    permitted temporal periodicity, and declared-record bounds.
14. Configuration, Nyquist, temporal-alias, memory and counter-range guards.
15. Record bounds at a large absolute time origin.
16. Rejection of normalization that cannot retain a meaningful representable
    factor, with a valid small-factor control case.
17. Gaussian tails at wavenumber scales `1e-200` and `1e200` against a closed-form
    expected value.
18. Proper complex Gaussian covariance, pseudo-covariance, mean and fourth
    moment over 4096 spaced seeds.
19. Physical-field covariance and the mean/variance of realized kinetic energy
    over another 4096 spaced seeds.

The full-lattice tensor check is independent of the implementation's
half-lattice cross-product factorization. Tests use approximately `2e-13`
scaled absolute tolerances for algebraic comparisons, `1e-12` normalized mode
transversality, and `5e-9` scaled central-difference tolerances. These are core
tests and do not establish a spectral/hp solver's discrete divergence accuracy.

Statistical acceptance bands use six analytically derived standard errors.
The Gaussian fourth-moment check distinguishes Gaussian amplitudes from a
fixed-magnitude random-phase substitute. The realized-energy variance uses the
quadratic Gaussian energy's moments, not a Gaussian assumption for energy.
Seeds are widely spaced: the frozen `seed+6*p+q` addressing rule would reuse
shifted raw-word windows for consecutive integer seeds, so treating those as
independent ensemble realizations would miscalibrate uncertainty.

Independent review also caught two issues fixed before the tested commit: an
absolute-time tolerance that could grow larger than its record, and a positive
normalization rounding to zero in diagnostic metadata. Both have explicit
regressions. The record tolerance is capped relative to its duration, and an
unrepresentable normalization fails rather than producing inconsistent metadata.

## Linux build and numerical evidence

The [machine-readable receipt](hit-inflow-milestone-2-validation.json) records
the exact tested revision, nine source/build-file digests, fixture-manifest
digest, installed-header digests, commands, settings and full-log digests.

The clean Linux checkout at `/media/jay/新加卷/nektar-v2026` was synchronized to
`a6790a4a505a78e363759383478599efe2bd837b` before compilation. Build and install
used `/home/jay/build/nektar-v2026-mpi` and
`/home/jay/.local/opt/nektar-v2026-mpi`. Logs and temporary files are under
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-2/a6790a4`.

Configuration: GCC 11.4.0, C++17, Release with `-O3 -DNDEBUG`, MPI enabled,
FFTW disabled, unit tests enabled. The core needs no FFTW backend. Compilation
used at most 16 workers; tests used one serial process, no MPI ranks, and
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`. There were no active compiler, solver,
or test-runner jobs at the allocation check. Commands ran sequentially and
remained within the 72-worker total limit.

After sourcing the installation's `env.sh` and setting those thread limits:

```bash
cmake -S '/media/jay/新加卷/nektar-v2026' \
    -B /home/jay/build/nektar-v2026-mpi
cmake --build /home/jay/build/nektar-v2026-mpi \
    --target SolverUtilsUnitTests --parallel 16
ctest --test-dir /home/jay/build/nektar-v2026-mpi \
    -R '^HITCoreUnitTests$' --output-on-failure --parallel 1 -V
```

Result: all 19 Boost cases passed; the focused CTest passed 1/1 in approximately
0.11 seconds including its test-runner overhead. A separate diagnostic capture
selected the three ensemble-related cases with `--log_level=message` and
`--report_level=detailed`; all 68 assertions passed. Its rounded output was:

| Quantity | Sample result | Theoretical target |
| --- | --- | --- |
| Complex Gaussian component `E[abs(xi0)^2]` | 1.01875 | 1 |
| Mean realized kinetic energy | 0.0739417 | 0.073926 |
| Variance of realized kinetic energy | 0.00026885 | 0.00026893 |

These passed the predeclared sampling-error bands; a finite sample is not
expected to give the target moments exactly. The reported values are rounded
diagnostics, while the tests use their full double precision.

`cmake --install /home/jay/build/nektar-v2026-mpi` succeeded. Both installed
public headers under `include/nektar++/SolverUtils/Turbulence/` match the tested
source bytes. No incompressible solver simulation or downstream experimental
comparison was run at this core milestone.

## Closure and next step

Milestone 2's completion criteria are met: both spectrum options, validation,
cutoff/taper, deterministic Gaussian coefficients, conjugate pairing, zero-mode
handling, represented-energy diagnostics, and a reference evaluator with
analytical derivatives are implemented and validated on Linux. Useful errors
cover malformed inputs, overflow, invalid normalization and resource limits.

Milestone 3 can now construct the full-3D inlet class using these immutable
modes. It must implement and test the all-rank configuration load and the
prescribed-velocity hook before pressure acceleration is formed. MPI
distribution, restart/checkpoint serialization, and solver time-accuracy claims
remain subject to the later integration milestones.
