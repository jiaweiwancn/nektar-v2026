# Neutral ABL RANS kernel and wall-layer tests

Milestone 5 adds `RANSWallModel=AWFEquilibrium` for CLS. Fixed weighting
isolates the dissipation correction; switched weighting uses seven independent
resolved indicators and explicitly supplied `RANSAWFOn`/`RANSAWFOff`
thresholds. The current configured CLS reference, cutoff admissibility,
original-form fallback and strict preflight are described in the
[Milestone 5 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-5.md).
Use `AWFEquilibrium` Robin tags for tangential velocity and k. The small
`AWFEquilibrium.xml` fixture is registered as `RANSAWFEquilibriumSmoke` in
the production and external builds. Run the complete Linux study with:

```bash
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m5.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/validation-m5"
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/plot_m5.py \
  "$RANS_RUN/validation-m5/validation.json" "$RANS_RUN/figures-m5"
```

The study audits actual native stage residuals, full fields, surface budgets,
and per-point switch/source/traction receipts against independent quadrature.
Its `--pilot` option cannot close the milestone. Cases are reused only when
the exact session and both solver/audit executable fingerprints match.
All numerical threads are one and cases run sequentially, using at most
four MPI ranks. Raw point CSVs and case receipts remain under the Linux run
directory; the compact validation manifest records their hashes and summaries.
No default switch thresholds imply terrain calibration. Full inlet/outlet,
restart and longer preservation studies remain Milestone 6.

Milestone 4 adds `TurbulenceModel=CLSHighRe` with primary dissipation
`eps_tilde`, the full cubic covariance and algebraic Yap source. The optional
`RANSWallModel=AWFOriginal` couples the original wall function through the
conservative CV-to-CG strip. Its [specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-4.md)
defines the physical frame, forcing, strip inertia, normal velocity, algebraic
dissipation, ownership, and conservation diagnostics. `CLS_MMS.xml` and
`AWFOriginal.xml` are short CTest fixtures. With the isolated full executable
built as below, run the sequential Linux acceptance study with

```bash
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m4.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/validation-m4"
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/plot_m4.py \
  "$RANS_RUN/validation-m4/validation.json" "$RANS_RUN/figures-m4"
```

It uses at most four MPI ranks, one case at a time, with numerical threads set
to one. The generator uses SymPy for independent manufactured forcing and
SciPy quadrature for the original-AWF oracle. Original-AWF wall data use Robin
tangential/k entries tagged `AWFOriginal`, Dirichlet normal velocity and
positive initial epsilon data, and the coupled pressure marker. Epsilon and
normal traces are subsequently determined within the nonlinear stage.
The original baseline intentionally retains the note's wall-source imbalance.
The equilibrium correction and switch are added by Milestone 5.

Milestone 3 provides an opt-in `RANSVelocityCorrectionScheme` for standard
k–epsilon URANS. Its [specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-3.md)
defines the implicit nonlinear stage, mixed pressure space, boundary fluxes,
pressure convention, and scope. Build the complete executable in isolation
by adding `-DRANS_BUILD_SOLVER=ON` to the CMake command below. It compiles the
production incompressible-solver source list against the installed native
libraries and preserves the shared installation used by other experiments.

Use `StandardRANS_MMS.xml` for a short manufactured smoke test, registered as
`RANSStandardTransportSmoke` in both the main and external builds. The full
Milestone 3 study is:

```bash
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m3.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/validation-m3"
```

Set `RANS_RUN` and TMPDIR as described below. The fixture generator needs SymPy;
the validator otherwise uses the Python standard library and the native
`RANSFieldAudit` executable. Temporal cases use at most three concurrent
one-rank jobs. That pool closes before sequential four-rank MPI cases start.
Affinity and library threads are limited before loading generator dependencies.
`--skip-long-abl` is a pilot option and cannot close the milestone. Completed
cases are reused only when their exact fixture and executable hashes match.

The accepted configuration uses full 2D/3D CG velocity and scalar fields, with
degree p−2 continuous pressure at the same quadrature points. All physical
pressure boundaries use the `RANSPressure` Neumann marker with prescribed normal
velocity. The lower standard wall uses tangential `RANSLogLaw` Neumann data,
zero normal Dirichlet velocity, zero k flux and algebraic epsilon Dirichlet
data at a positive standoff height. See the supplied fixture and generator for
complete XML, including explicit quadrature settings. Saved `p` means reduced
pressure, `pi = p_mechanical/rho + 2k/3`.

Milestone 1 selects a control-volume wall strip on `0 < z < H`, coupled to
native continuous spectral elements on `H <= z <= Z`. The mathematical
contract and completion evidence are in
[the Milestone 1 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-1.md).

`RANSWallLayerProbe` checks the separate wall rows, dissipation constraint,
normal mass balance, interface signs, native outer weak residuals, h/p and
quadrature convergence, inconsistent-boundary controls, and a conservative
backward-Euler step with frozen reference viscosity. The outer residual now
calls the general CLS kernel at the native simple-shear gradients, including
its variable coefficient, full stress, and algebraic Yap term. This probe
isolates the wall strip; native coupled transport is covered by Milestone 4.

Milestone 2 adds `RANSKernelProbe` and an independent Python oracle. They test
standard k–epsilon and CLS tensors/sources, configured reference roots, original
AWF branches and forcing, corrected dissipation budgets and profile consistency,
admissibility, baseline recovery, and independent switch sensors. Full contracts
and measured results are in
[the Milestone 2 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-2.md).

The fixed reference is the note's rough-wall example with `P=1`, `H=2`, and an
additional matching sample at `N=4` metres. The CG intervals are `[2,4]`,
`[4,8]`, `[8,16]`, and `[16,30]` metres before uniform subdivision of each
interval. Geometry and sample locations do not change with polynomial order.
The N sample belongs to the element on its right, exercising communication
from rank 1 to the wall-owning rank 0 in the four-rank test.

The executables are registered as `RANSWallLayerPrototype` and
`RANSKernelUnitTests` in the main build when `NEKTAR_BUILD_TESTS=ON`.
An external build compiles the same production kernel sources into a static
library and links the installed native Nektar++ libraries, preserving the
shared build and installed binaries:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export RANS_RUN=/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/m2-reproduction
mkdir -p "$RANS_RUN/tmp"
export TMPDIR="$RANS_RUN/tmp"
taskset -c 112-127 cmake -S solvers/IncNavierStokesSolver/Tests/RANS -B "$RANS_RUN/build" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DCMAKE_BUILD_TYPE=Release
taskset -c 112-127 cmake --build "$RANS_RUN/build" --parallel 4
taskset -c 112-127 ctest --test-dir "$RANS_RUN/build" -R '^RANS(WallLayerPrototype|KernelUnitTests)$' \
  --output-on-failure --parallel 1
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m2.py \
  --build "$RANS_RUN/build" \
  --output "$RANS_RUN/validation" --mpi-ranks 4
```

Run from the synchronized Linux source root. Choose a fresh experiment directory;
the commands above use a separate reproduction directory. CPU IDs 112–127 are
the recorded allocation on this server; use up to 16 allowed IDs if its affinity
configuration changes. The project limit is **16 CPU workers in
total**, including concurrent build and test work. The validator serializes its
one-process and four-rank runs, sets library thread counts to one, and limits
its inherited Linux CPU affinity to at most 16 logical CPUs. Build and CTest
allocations must also share the 16-worker budget. The repository-wide limit
across other server projects still applies.

Outputs are the kernel report, serial/MPI column reports, their logs, and a
validation manifest with source hashes, commands, CPU allocation, thread
settings, and comparisons. Milestone 2's oracle requires NumPy and SciPy
(tested with 1.21.5 and 1.8.0). It sets BLAS/thread environment variables and
CPU affinity before importing those packages. The kernel probe itself has no
NumPy/SciPy dependency.

For the original standard-library-only column audit, use `validate_m1.py`
with `--probe "$RANS_RUN/build/RANSWallLayerProbe"`, `--output` and
`--mpi-ranks 4`. It now exercises the general CLS kernel through the column.
Keep historical Milestone 1 reports unchanged when running this regression.

A successful kernel/column result establishes the specified local model
contracts and column convergence. Milestone 3 adds standard-model unsteady
transport, Milestone 4 adds CLS and the original AWF, and Milestone 5 adds the
corrected-AWF coupling. The full analytical inlet/outlet preservation study
remains a later milestone.
