# H14 fixed-height mapped backend

`MappedBlockV1` is an opt-in serial route for the qualified horizontal,
zero-hill `Plane` / `CompositeTunnelV1` / `OpenTunnelV1` WJ Budget case.
The default `Dense` route and its 3000-state guard remain. Its older 4331-state
exception still excludes mapped cases. H14 uses only the stored H13 6 mm case;
the hill partner remains geometry evidence, without governing evolution.

The new route stores field operators sparsely, permutes their actual graphs,
and factors consistent mass with LAPACK band Cholesky and nonsymmetric blocks
with band LU. The pressure-only nonsymmetric Schur matrix is bounded at 1024
rows. Columns of `B G` are streamed and discarded, then the correction and true
Schur residual are independently applied. Here `B=(alpha M_free)^-1`; it is
neither a lumped mass nor the momentum inverse. Bound whitening and recovery
use the same permuted consistent-mass factor and retained-pool certificates.

Selection and capability rejection occur before patch mass/reference allocation.
Targets above 3000 states also need `RANSWallMappedResourceDesign`, the SHA-256
of their admitted immutable resource design. The finite state limit is 24688,
with a 512 MiB per-factor payload cap checked before allocation. This is a
transitional serial backend, not a scalable production solver.

## Build and focused checks

On Linux, use the environment and qualified prebuilt SpatialDomains library
recorded in the H14 evidence. Build both the governor and strict reader from
the same source; do not replace a shared installation. Example from this run:

```bash
R=/home/jay/build/wj-h14-mapped-backend-20261113
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cmake -S "$R/source/solvers/IncNavierStokesSolver/Tests/RANS" -B "$R/build" \
  -DCMAKE_BUILD_TYPE=Release -DRANS_BUILD_SOLVER=ON \
  '-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math' \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_PREBUILT_SPATIAL_LIBRARY=/home/jay/build/wj-h10-flat-inlet-20261113/runtime/geometry/libSpatialDomains.so.5.9.0
cmake --build "$R/build" --parallel 8 --target \
  IncNavierStokesSolver RANSWallResolutionProbe WJMappedAlgebraChecks RANSKBoundQPChecks
```

`WJMappedAlgebraChecks` covers restrictions/permutations, nonsymmetric and
transpose solves, streamed pressure algebra, retained epsilon columns, and
activated/near-dependent bounds. `RANSKBoundQPChecks` retains the existing QP
controls. The executed commands/logs, rather than registration count, determine
which tests passed. `test_wj_h14_execution.py` has twelve executed methods for
the execution diff, immutable reference, owner index/cache, resource halt,
completed-state selection and bounded native-packet consumer.

Use qualified runtime wrappers for physical probes. H14 copied the original
private H10 loader/libc closure and verified actual native-exec maps against its
sealed dependency inventory. Directly launching a new ELF with host libc is
not the qualified runtime.

## Case and observation entrypoints

`wj_h14_execution.clone` applies execution-only changes to the manifest-bound
H13 XML; it does not generate a mesh, fit an intercept or invert roughness.
`admit_wj_h14.py --mode prepare` requires complete small qualification,
preflights actual spaces/allocations, and binds the new resource design.
The pilot timestep is frozen at .0025. The stored H13 physical Functions,
normal partition, matching height and observation roles are preserved.

`run_wj_h14.py` is scoped to this root and its persistent ledger. Reservations
are write-once, and scientific halts are monotone. `--phase target_cp0` and
`--phase pilot` cannot be repeated. This published run is closed; these commands
are reproduction documentation, not authorization to reuse its allocation.

`study_wj_h14.py --mode assessment --case target/pilot1 --new-only` assesses
only genuinely new completed states and uses the actual CP0 for drift. A restored
alpha-zero entry is not a completed callback. The owner inventory includes every
legitimate face owner and physical Gauss8 support. Geometry/basis caches contain
immutable coordinates only; current K, cutoff, stress and reactions are refreshed.

For strict reading, `--mode view` writes an `open_checks` intent with its genuine
producer/checkpoint/history. Run the qualified reader once, then use
`--mode view-finish` on the resulting bytes. `weak_compact=true` retains native
projected `1`, `x/5`, `interface_hat` signed components/test hashes. This performs
no new independent quadrature campaign. `compact_open_packet` decodes the
summary without inflating millions of retained scalar observation rows; the
complete native file remains immutable and hash-bound.

Mapped k snapshots use schema 2 COO records `(int32 row, int32 column, binary64
value)`. Their row maps, shape, encoding, byte count and hash are explicit;
legacy schema-1 dense readers must refuse rather than reinterpret them.
Legacy dense arithmetic/matrix observers and verification modes are rejected
before mapped allocation. They are not silently truncated.

The opt-in physical lookup prioritizes exact element containment before a
roundoff-only neighbor. An exact-face tie still uses the established deterministic
owner rule. H14 retains the original failed one-sided pressure-action witness
and the passing one-ULP normal/horizontal owner regression.

## Published outcome and next admission

One fresh 24688-state CP0 and one accepted, qualified full-domain step to .0025
were executed. The remaining three potential pilot steps did not fit the same
600-second pilot reservation with 25% headroom and were not launched. No duration
or hill experiment followed. See the H14 evidence for exact identities, factors,
costs, failures, strict residuals and closed accounting.

Before a separately scoped preservation experiment, qualify persistent immutable
row/test caching and a bounded warm-block cost/monitoring policy. Re-admit the
unchanged H13 clocks, time/error budgets and required observation window using
measured costs. Do not interpret this pilot as temporal convergence, equilibrium
preservation, physical height selection or experimental validation.
