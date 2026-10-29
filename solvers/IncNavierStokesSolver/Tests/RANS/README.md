# Neutral ABL wall layer prototype

Milestone 1 selects a control-volume wall strip on `0 < z < H`, coupled to
native continuous spectral elements on `H <= z <= Z`. The mathematical
contract and completion evidence are in
[the Milestone 1 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-1.md).

`RANSWallLayerProbe` checks the separate wall rows, dissipation constraint,
normal mass balance, interface signs, native outer weak residuals, h/p and
quadrature convergence, inconsistent-boundary controls, and a conservative
backward-Euler step with frozen reference viscosity. The outer residual uses
the simple-shear reduction of CLS, including its variable coefficient and
algebraic Yap term. The general cubic closure and the complete URANS timestep
are subsequent milestones.

The fixed reference is the note's rough-wall example with `P=1`, `H=2`, and an
additional matching sample at `N=4` metres. The CG intervals are `[2,4]`,
`[4,8]`, `[8,16]`, and `[16,30]` metres before uniform subdivision of each
interval. Geometry and sample locations do not change with polynomial order.
The N sample belongs to the element on its right, exercising communication
from rank 1 to the wall-owning rank 0 in the four-rank test.

The executable is registered as `RANSWallLayerPrototype` in the main build when
`NEKTAR_BUILD_TESTS=ON`. An external build can also link the already installed
Nektar++ libraries, preserving the shared build and installed binaries:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S solvers/IncNavierStokesSolver/Tests/RANS -B "$RANS_RUN/build" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$RANS_RUN/build" --parallel 4
ctest --test-dir "$RANS_RUN/build" -R '^RANSWallLayerPrototype$' \
  --output-on-failure --parallel 1
python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m1.py \
  --probe "$RANS_RUN/build/RANSWallLayerProbe" \
  --output "$RANS_RUN/validation" --mpi-ranks 4
```

Set `RANS_RUN` to a dedicated experiment directory under
`/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf`, create it first, and
set TMPDIR to a directory there. The project limit is **16 CPU workers in
total**, including concurrent build and test work. The validator serializes its
one-process and four-rank runs, sets library thread counts to one, and limits
its inherited Linux CPU affinity to at most 16 logical CPUs. Build and CTest
allocations must also share the 16-worker budget. The repository-wide limit
across other server projects still applies.

Outputs are serial/MPI JSON reports, their logs, and a validation manifest with
source hashes, commands, CPU allocation, thread settings, and comparisons.
No external Python packages are needed. A successful result establishes the
specified prototype contracts and convergence, not exact finite-order
preservation of the full nonlinear analytical inlet in IncNavierStokesSolver.
