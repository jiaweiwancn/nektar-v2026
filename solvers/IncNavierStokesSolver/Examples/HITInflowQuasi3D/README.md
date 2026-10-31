# HIT inflow with a homogeneous Fourier direction

This complete session uses four xy quadrilaterals, degree 7 and eight full real
Fourier slots in z (`HomModesZ=8`). It evolves all three velocity components
using `Homogeneous1D`; it is not a two-component 2D turbulence model. The domain
is 1 by 2*pi by 2*pi m, periodic in y and Fourier-periodic in z. The x-minimum
inlet prescribes HIT velocities and pressure `H`; the outlet uses `HOutflow`.

`shell-spectrum.dat` is an explicitly labelled manufactured shell spectrum,
with linear interpolation, zero extensions and `TargetRMS=0.05 m/s`. It has
no experimental values. The lattice (4,4,4), seed 20261028 and 1.55/1.5 rad/m
cutoff/taper retain a very small mode set. Ubar=1 m/s, nu=0.01 m2/s,
dt=0.001 s and 20 steps make this a configuration demonstration. It starts from
the mean flow. No statistical, DNS or experimental validation follows from this
short run. Spectral/hp dealiasing is enabled; homogeneous dealiasing and SVV
are off in the default matrix-transform session.

## Run in a fresh directory

The session omits `USEFFT`, selecting native matrix transforms. It therefore
does not require an FFTW build. This project-server command uses two MPI ranks:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
hit_src=/media/jay/新加卷/nektar-v2026
hit_run=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/examples/h1d-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$hit_run/tmp"
cp "$hit_src/solvers/IncNavierStokesSolver/Examples/HITInflowQuasi3D/HITInflowQuasi3D.xml" "$hit_run/"
cp "$hit_src/solvers/IncNavierStokesSolver/Examples/HITInflowQuasi3D/shell-spectrum.dat" "$hit_run/"
cd "$hit_run"
export TMPDIR="$hit_run/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mpirun --bind-to none -n 2 IncNavierStokesSolver --npz 2 HITInflowQuasi3D.xml > solver.log 2>&1
```

For serial execution in another fresh directory, use
`IncNavierStokesSolver HITInflowQuasi3D.xml`. An MPI-enabled build is needed for
the two-rank command. `HomModesZ/npz` must be even; the synthesis count `NZ`
is separate from `HomModesZ`. Keep total concurrent Linux workers at or below
64, with one OpenMP/BLAS thread per rank.

To use FFTW, add `<I PROPERTY="USEFFT" VALUE="FFTW" />` in `SOLVERINFO` in a
fresh copy and use an installation built with `NEKTAR_USE_FFTW=ON`. Matrix and
FFTW paths use the same mode/time realization; their native transform roundoff
can differ. Native homogeneous dealiasing (`DEALIASING=True`) additionally
requires `USEFFT=FFTW`; enabling it changes the numerical configuration and
should be checked separately.

The solver should finish 20 steps at t=0.02 s. Outputs include
`HITInflowQuasi3D.fld`, checkpoints, an `.opt` file and `history.his` with 12
physical probes and 21 frames including t=0. With optional test tools installed,
a serial diagnostic can read the completed parallel field directory:

```bash
HITFlowDiagnostics HITInflowQuasi3D.xml HITInflowQuasi3D.fld > diagnostics.log 2>&1
```

It writes `HITInflowQuasi3D.flow.json`. Do not treat its finite numerical
divergence, the generic exact-solution norms, or 20 steps of probe data as a
qualified turbulence-decay result.

## Measured spectra and comparison reproduction

The input header declares `KIND=Shell3D`, angular k in rad/m and E in m3/s2.
The loader rejects a `Longitudinal1D` spectrum as inlet data. The
[HIT data guide](../../Tests/HITInflow/README.md) gives exact SI preparation of
the upstream CBC Table-3 shell spectrum and keeps downstream measurements
separate. Replacing this toy table with measurements requires new resolution,
cutoff, normalization and recurrence checks; it is not an equivalent example.

The [full-3D example](../HITInflow3D/README.md) demonstrates the von Karman input.
The [pilot guide](../../Tests/HITInflow/pilots/README.md) documents production
settings and regeneration of the existing M7 comparisons. M7's full-band
agreement and 2% precision were not established; narrow-band validation remains
deferred.

An initial legacy `_0.chk` at t=0 can precede HIT metadata setup and cannot be
used as a verified HIT restart when identity metadata is absent. The evolved
checkpoints at t=0.01/0.02 s contain the identity; retain all four fields and
check the stored time before restarting.
