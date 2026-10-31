# HIT inflow examples, tests and experimental fixtures

This directory contains the Milestone 1 configuration examples and the principal
5.08 cm grid dataset from Comte-Bellot and Corrsin (1971),
[Journal of Fluid Mechanics 48, 273-337](https://doi.org/10.1017/S0022112071001599).
The implemented `HITInflow` condition supports full 3D and Homogeneous1D
three-component velocity correction. The `hit-inflow-*.xml` files here are
configuration fragments. Complete runnable examples are in
[HITInflow3D](../../Examples/HITInflow3D/README.md) and
[HITInflowQuasi3D](../../Examples/HITInflowQuasi3D/README.md).
The [user guide](../../../../docs/user-guide/solvers/incompressible-ns.tex)
contains the version-1 configuration reference (`sec:hit-inflow`).

The normative interface and integration contracts are in
[the Milestone 1 specification](../../../../docs/plans/hit-inflow-milestone-1.md).
The original project sequence is in
[the implementation plan](../../../../docs/plans/hit-inflow-implementation-plan.md).

## Source and units

`cbc1971-source.json` preserves the checked source numbers as decimal strings
in centimetre units, including their printed precision and missing entries.
The complete Table 2(a), Table 3, and principal Table 4 rows were checked
visually against the local original PDF on 2026-10-28. Their printed pages are
296, 298, and 299 (PDF pages 24, 26, and 27). The source record contains that
PDF's SHA-256 digest. The source PDFs are supplied separately under
`docs/literature/` and are not required to regenerate the numeric fixtures.

The measurement-region speed is 12.7 m/s, from the footnote on printed page
289 (PDF page 17). The approach speed is 10 m/s. Stations use
`theta = U0*t/M`, with travel time through the upstream contraction, rather
than simply `x/M`. Taking station 42 as the computational reference gives
station offsets of 3.612896 m and 8.322564 m at stations 98 and 171.

Wavenumbers are angular wavenumbers. Multiply source `cm^-1` by 100 to obtain
`rad/m`. Multiply either spectrum in `cm^3/s^2` by `1e-6` to obtain `m^3/s^2`.
No extra factor of `2*pi` belongs in this conversion.

## Numeric files

| Files | Meaning | Role |
| --- | --- | --- |
| `cbc1971_theta42_spectrum_3d.dat` | Table 3 shell spectrum, 19 points | Experimental inlet input |
| `cbc1971_theta98_spectrum_3d.dat` | Table 3 shell spectrum, 19 points | Downstream validation |
| `cbc1971_theta171_spectrum_3d.dat` | Table 3 shell spectrum, 18 points | Downstream validation |
| `cbc1971_theta*_longitudinal_1d.dat` | Table 2(a) one-sided longitudinal spectra, 22 points each | Virtual-probe validation |
| `cbc1971_stations.dat` | Table 4 converted to SI, with relative times and distances | Gross statistical targets |

The first comment in each spectrum declares its kind and units. Both kinds have
two numeric columns, but their meanings differ. The inlet loader accepts `KIND=Shell3D` and rejects `KIND=Longitudinal1D`. Shell energy integrates to
one half of total velocity variance; the one-sided longitudinal spectrum
integrates to streamwise velocity variance.

Missing source entries remain missing: individual station files omit them.
There is no smoothing, extrapolation, taper, or normalization in these files.
The small nonmonotonic sequence in Table 3 at station 98 is printed in the
original and is preserved. The source spectrum has finite coverage, and matching
its measured part does not establish complete energy or dissipation integrals.
The runtime contract makes extensions explicit.

The Table 4 consistent initial viscosity is about `1.5e-5 m^2/s`. Figure 9's
caption instead uses `2*nu=0.28 cm^2/s`; this source discrepancy remains recorded.
`lambda` retains the paper's Taylor transverse microscale convention. The last
station-file column is the dimensionless product `(lambda/Lg)*Re_lambda`.

## Reproduction and checks

From the repository root, using only Python's standard library:

```bash
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/prepare_cbc1971.py
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/prepare_cbc1971.py --check
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/check_contract.py
```

The preparation command rewrites only the seven generated `.dat` files and
`cbc1971-manifest.json`. It uses exact decimal scaling and writes a deterministic
manifest with source, generator, and output SHA-256 digests. `--check` performs
no writes and fails on changed data or digests. `--output-dir` permits comparison
against a separate directory. To verify a particular supplied source copy too:

```bash
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/prepare_cbc1971.py --check \
  --source-pdf 'docs/literature/Simple Eulerian time correlation of full-and narrow-band velocity signals in grid-generated isotropic turbulence.pdf'
```

`check_contract.py` checks these examples and their fixture-level constraints.
It is not the future C++ parser and does not certify a solver run. The numerical
generator, pressure convergence, MPI sampling, and physical comparisons have
their own later milestone checks.


## Solver checks and completed comparison

Small full-3D, Homogeneous1D, MPI, transform-backend and restart checks are
documented in [README-m3.md](README-m3.md), [README-m4.md](README-m4.md) and
[README-m5.md](README-m5.md). The solver/test build needs
`NEKTAR_BUILD_TESTS=ON`, which also builds the native diagnostic tools.
`NEKTAR_BUILD_UNIT_TESTS=ON` additionally enables the core/boundary/restart/SVV
unit groups. The project installation has MPI and FFTW enabled.
Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`, and keep all
concurrent work within the current 64-worker Linux limit. Use a fresh directory
under `/media/jay/新加卷/nektar-v2026-runs/hit-inflow/` for generated output.

The [pilot guide](pilots/README.md) describes the frozen production recipe,
source-data preparation, analysis requirements and the sequential
`pilots/reproduce_m7.py` command that regenerates all four existing analyses
and five figure pairs without a new CFD run. The
[M7 report](../../../../docs/plans/hit-inflow-milestone-7.md) and
[receipt](../../../../docs/plans/hit-inflow-milestone-7-validation.json) retain
the tested native revision, MPI/thread settings, input hashes and data locations.
Both spectrum models decay too slowly and retain too much correlation; all
32 sampling screens fail the 2% criterion. Full-band physical validation is not
established, and narrow-band validation is deferred.

## Reproduce the configured small regressions

The following exports the configured CTest commands into separate fresh working
folders rather than writing generated fields beside source inputs. It runs one
test at a time. With MPI, FFTW and optional unit tests enabled, the selection is
14 solver regressions plus four native unit groups (48 cases). The largest
solver regression uses eight MPI ranks. Inspect current registrations with
`ctest --test-dir "$hit_build" -N -R HIT` before execution.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
hit_build=/home/jay/build/nektar-v2026-mpi
hit_check=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/regression-reproduction-$(date -u +%Y%m%dT%H%M%SZ)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
taskset -c 0-7 python3 - "$hit_build" "$hit_check" <<'PY'
import json, pathlib, subprocess, sys
build, output = sys.argv[1], pathlib.Path(sys.argv[2])
tests = json.loads(subprocess.check_output([
    'ctest', '--test-dir', build, '--show-only=json-v1', '-R', 'HIT'], text=True))['tests']
output.mkdir(parents=True, exist_ok=False)
for test in tests:
    properties = {p['name']: p['value'] for p in test.get('properties', [])}
    if int(properties.get('PROCESSORS', 1)) > 8:
        raise RuntimeError('test exceeds this eight-worker allocation')
    work = output / test['name']
    work.mkdir()
    with (work / 'test.log').open('x') as log:
        subprocess.run(test['command'], cwd=work, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    print('PASS', test['name'], flush=True)
PY
```

Reserve those eight CPUs from the combined 64-worker budget before running.
Native MPI helper threads inherit the selected CPU set; OpenMP/BLAS remain
one per rank. The exact registered commands and individual logs are retained in
the M8 closure receipt/run directory. Build options control which optional tests
are listed; the documented counts refer to the tested project installation.
