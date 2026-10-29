# Neutral ABL RANS kernel and wall-layer tests

Milestone 1 selects a control-volume wall strip on `0 < z < H`, coupled to
native continuous spectral elements on `H <= z <= Z`. The mathematical
contract and completion evidence are in
[the Milestone 1 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-1.md).

`RANSWallLayerProbe` checks the separate wall rows, dissipation constraint,
normal mass balance, interface signs, native outer weak residuals, h/p and
quadrature convergence, inconsistent-boundary controls, and a conservative
backward-Euler step with frozen reference viscosity. The outer residual now
calls the general CLS kernel at the native simple-shear gradients, including
its variable coefficient, full stress, and algebraic Yap term. The complete
URANS timestep remains a subsequent milestone.

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

A successful result establishes the specified local model contracts and column
convergence. Full nonlinear URANS integration and analytical inlet preservation
in IncNavierStokesSolver remain later milestones.
