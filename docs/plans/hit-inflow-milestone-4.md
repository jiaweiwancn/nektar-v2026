# HIT inflow Milestone 4 closure: Homogeneous1D, MPI and restart

Milestone 4 is **CLOSED on 2026-10-29**. The common HIT generator now supplies
full-3D and multiple-mode Homogeneous1D inlets through native constrained
boundary transforms. Linux validation passes for matrix and FFTW transforms,
mesh and Fourier MPI partitions, cross-representation refinement, and verified
restart with preserved inlet phase.

This closes Milestone 4 of the [implementation plan](hit-inflow-implementation-plan.md).
It does not close the downstream-decay experiment or the production-resolution
study. All numerical conclusions below use the authoritative Linux environment.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | implementation and numerical verification |
| Origin date | 2026-10-29 |
| Version label | hit_inflow_homogeneous_mpi_v1 |
| Verification status | VERIFIED for the stated Linux cases and acceptance bounds |
| Final implementation revision | `bbc49bd91d1e3b7e282a4dd9edae7de3de2258e5` |
| Branch | `codex/hit-inflow-generator` |
| Numerical receipt | [Exact revisions, commands, hashes and results](hit-inflow-milestone-4-validation.json) |

The final build, 35 unit cases and all restart checks use the final revision.
The baseline MPI/backend batch uses runtime `8bd06ea` with final fixture
settings from `0049141`; polynomial/Fourier refinement uses `2b33286`.
The ten selected regressions use `8bd06ea`. Subsequent runtime changes only
verify source checkpoint metadata and prohibit a restart-time phase reset.
They do not change synthesis, sampling, transforms, mesh serialization, flow
equations, integration timing, or ordinary expression-based startup. The receipt
records full revision IDs rather than treating every run as one checkout state.

## Implemented behavior

- `HITInflowConfig` accepts full 3D and full-Fourier Homogeneous1D, validates
  the homogeneous period and strict cutoff below the solver's z Nyquist limit,
  and rejects single/half modes and homogeneous strips.
- `HITInflowGeometry` checks four xy box sides in Homogeneous1D and six sides
  in full 3D. It counts the extruded boundary area once, on the owner of global
  plane zero, and handles ranks without inlet faces through shared collectives.
- The existing boundary sampler uses global physical-plane coordinates,
  fills all three components with mean plus fluctuations, sets physical-space
  state and calls `FwdTransBndConstrained` once. That native operation already
  includes the homogeneous forward transform; there is no extra transform.
- Synthesis `NZ` defines the random lattice independently of solver
  `HomModesZ`. Increasing solver planes can refine one fixed realization.
  Both synthesis and solver Nyquist restrictions still apply. This clarification
  replaces the original `NZ=HomModesZ` restriction in the M1 specification.
- `HITInflowRestart` records 19 required identity entries: schema, algorithm,
  seed, periods/counts, origin, mean, time origin, spectrum definition and source
  digest, cutoff/taper, normalization/RMS/scale, numerical environment,
  normalized table, coefficient digest and overall identity digest.
- Restart verification reads source metadata into an empty map. It requires
  stored time, matching identity and a valid complete job interval before the
  first inlet trace. A native `--set-start-time` override may only restate the
  stored time within floating-point tolerance. Duration, repeat policy, output
  schedule, timestep and MPI layout remain independent validity checks.
- The preferred `Checkpoint` filter now preserves equation metadata and writes
  actual time at full floating-point precision, alongside its checkpoint index.
  Legacy checkpoints and final fields already carry equation metadata.

Two native tools support the checks: `HITBoundaryProbe` samples/inverts the
actual constrained trace on every rank, including repeated and backward times;
`HITFieldCompare` imports complete field directories and reports full-precision
physical norms/differences. The latter verifies representation, uniform expansion
modes, period and complete element/Fourier-plane coverage before comparison.

## Validation inputs and fixed criteria

The committed [M4 fixtures](../../solvers/IncNavierStokesSolver/Tests/HITInflow/README-m4.md)
provide 13 small sessions and a manufactured shell-spectrum table. The standard
box is x in `[0,1]`, y and z in `[0,2*pi]`, with four xy quads or four extruded
hexes, eight polynomial modes, eight Fourier planes, `dt=0.001`, and 20 steps
through `t=0.02`. Synthesis uses counts `(4,4,4)`, periods `(2*pi,2*pi,2*pi)`,
seed 20261028, mean speed 1, target RMS 0.05, cutoff 1.55 and taper start 1.5.
The von Karman case retains nine independent modes; the tabulated case retains
six oblique modes. These small bands test implementation behavior.

The numeric bounds in `m4-acceptance.json` are byte-identical to the pre-run
version at `6f8040b`. No acceptance bound was relaxed. Stochastic fixtures use
`IterativeStaticCond` with tolerance `1e-14`; exact uniform-flow fixtures use
`DirectStaticCond`. The acceptance file's note mentioning `1e-12` records the
initial trial tolerance, which proved insufficient for the full-3D MPI bound.

### Inlet and transform checks

All 18 probe/solver case pairs pass. Every probe's raw samples match direct
synthesis exactly, and repeated evaluation reproduces coefficients identically.
The largest serial projected error is below `1.9e-4` in full 3D and `7.6e-7`
in Homogeneous1D. Both spectra retain nonzero z coefficient energy in all three
components. The mean-only maximum L2/Linf error over velocity and pressure is
below `1e-10` in both representations.

For matching MPI/backend grids, raw sample differences are zero and identity
and coefficient digests agree. The table reports maximum projected inlet
component difference and evolved physical Linf differences from the native
field importer. Comparisons use the corresponding serial matrix case.

| Case | Projected inlet difference | Velocity Linf difference | Pressure Linf difference |
| --- | --- | --- | --- |
| H1D, mesh 4 | 0.000e+00 | 3.366e-12 | 1.392e-11 |
| H1D, Fourier 4 | 0.000e+00 | 4.022e-13 | 1.181e-11 |
| H1D, mixed 4 | 0.000e+00 | 3.369e-12 | 8.907e-12 |
| H1D, mixed 8 | 0.000e+00 | 3.369e-12 | 8.576e-12 |
| Full 3D, mesh 4 | 0.000e+00 | 1.642e-10 | 7.599e-10 |
| H1D table, mixed 4 | 0.000e+00 | 9.262e-13 | 1.712e-11 |
| H1D table, mixed 8 | 0.000e+00 | 9.182e-13 | 1.071e-11 |
| H1D FFTW, serial | 5.551e-16 | 2.309e-13 | 8.639e-12 |
| H1D FFTW, mixed 4 | 5.551e-16 | 3.374e-12 | 5.536e-12 |

The bounds are `5e-12` for matching-grid inlet samples, `2e-9` for evolved
velocity Linf, and `2e-8` for pressure Linf. In the mesh-four cases, ranks 0 and
3 own no inlet; in mixed-eight cases, ranks 0, 3, 4 and 7 own no inlet. Both
von Karman and tabulated Homogeneous1D cases pass with this ownership pattern.
Every written partition vertex was checked for exact double round trip against
the source mesh.

### Full-3D versus Homogeneous1D refinement

At common coordinates, raw inlet differences are at most `4.45e-16`, with the
same immutable identity. The von Karman projected-trace difference decreases
from `1.13388318e-4` at eight polynomial modes to `3.12565144e-6` at ten modes:
ratio `0.027566`, below the predeclared `0.1`. Both individual bounds,
`2e-3` and `5e-5`, pass. The tabulated eight-mode comparison also passes.

Full-precision norms are computed from native field data, separately from the
rounded solver error summary. The largest relative difference among velocity
component L2 norms decreases from `1.45059047e-4` to `2.89092857e-6`, a ratio
of `0.0199293`, below `0.5`; both are below `0.01`. For the tabulated eight-mode
case, this difference is `1.21655436e-4`. These are transient global flow norms
used to assess numerical agreement, not ensemble HIT statistics.

Increasing Fourier planes from 8 to 16 preserves the sampled realization and
projected inlet within `5e-12`. The maximum relative velocity-norm difference
is `3.27e-15` for this smooth small-band case; the pressure L2 norm difference
is `1.10e-11`. This does not establish sufficient Fourier resolution for a
production turbulent cascade.

### Restart and checkpoint checks

Both representations run continuously to `t=0.02` and as two segments split at
`t=0.01`. Stored timestamps are `0.010000000000000002`. With native IMEX history
reinitialized, the evolved restart difference is bounded rather than bitwise.
Halving dt from `0.001` to `0.0005` gives:

| Representation | Coarse maximum velocity Linf difference | Fine difference | Fine/coarse |
| --- | --- | --- | --- |
| Homogeneous1D | 1.15300949e-05 | 2.73235567e-06 | 0.236976 |
| Full 3D | 1.14858416e-05 | 2.72301139e-06 | 0.237075 |

All velocity differences are below `2e-4`, all pressure differences below
`2e-3`, and the refinement ratios are below `0.75`. The Checkpoint-filter
restart gives the same result as the legacy checkpoint path. A four-rank mesh
checkpoint restarted with mixed mesh/Fourier partitioning differs from the
serial restart by at most `2.85e-12` in velocity and `7.40e-11` in pressure.
An actual tabulated-spectrum serial checkpoint restarted on four mixed ranks
also passes, including the normalized-table XML round trip.

The guard suite has 23 passing outcomes, including 19 expected rejections.
It checks missing spectrum input; incompatible period, Nyquist limit, half
mode or strip; wrong field modes/period; missing Fourier-plane data; and an
otherwise valid single partition supplied instead of the full field directory.
Restart checks reject changed seed/environment, missing identity or time,
malformed time, phase-resetting start-time override, job overrun, mixed
file/expression initialization, different field files, and changed spectrum bytes
with unchanged numeric values. A valid duration/repeat-policy change passes.
MPI rejection cases exit with code 1 before stepping, within their timeouts.

### Unit and regression checks

All three focused CTest entries pass: 19 `HITCore`, 7 `HITBoundary`, and 9
`HITRestart` Boost cases, 35 total. The new source-metadata unit demonstrates
that equation defaults cannot supply missing stored time or identity, rejects
malformed/overridden time, and accepts a valid source map.

These ten solver regressions pass their unchanged metrics:

- `HITInflowMeanOnly3D`
- `HITInflowManufactured3D`
- `HITInflowHomogeneous1DMean`
- `ChanFlow_3DH1D_bcsfromfiles`
- `ChanFlow_3DH1D_FFT`
- `KovaFlow_3DH1D_MVM_FFTW_Consistency`
- `MovingRefFrame_translation`
- `MovingRefFrame_Rot_SimpleDomain`
- `ChanFlow_3DH1D_Parallel_mode1`
- `ChanFlow_3DH1D_Parallel_mode2`

## Diagnosed failures and corrections

The initial full-3D uniform-flow run at iterative tolerance `1e-12` exceeded
the `1e-10` pressure bound. Direct solves isolate the exact DC limit and pass.
The initial full-3D serial/MPI velocity difference was `5.24e-8`; tightening
the stochastic solve to `1e-14` reduces it to `1.64e-10`, with unchanged bounds.

Native XML partition writing originally stored vertices with eight digits
after the decimal in scientific notation: pi became `3.14159265`, moving inlet
coordinates by up to `3.6e-9`. The vertex writer now uses `max_digits10-1`
digits after the scientific decimal, retaining all significant double digits.
The MPI grid/sample checks and explicit vertex round trips pass after this fix.
Curved-edge and curved-face writers are unchanged.

The field validator initially compared the native single UNIORDER tuple with
an in-memory per-element tuple list, and required full precision from a native
six-significant-digit period header. It now accepts those exact native forms;
negative mode, period and complete-coverage checks still reject bad inputs.

A negative checkpoint with missing Time initially advanced from zero because
native metadata import merged into an equation map prepopulated with Time=0.
Fresh source-map verification fixes that defect. The subsequent time-override
guard enforces the existing no-phase-reset contract. Build/setup logs, failing
coordinate/tolerance trials and command history are retained separately from
accepted results in the Linux run storage and receipt.

## Linux environment and reproduction

Source: `/media/jay/新加卷/nektar-v2026`; build:
`/home/jay/build/nektar-v2026-mpi`; installation:
`/home/jay/.local/opt/nektar-v2026-mpi`. The run root is
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-4/6f8040b`.
The directory label identifies the initial candidate, not every tested revision.
Large/generated data stay outside the checkout. The receipt points to the
969-file raw-artifact manifest and records its SHA256.

This Release build has MPI and FFTW enabled and HDF5 disabled. FFTW development
headers were extracted from an Ubuntu package into isolated experiment storage;
the existing runtime library is used. Package, header and runtime-library
versions/hashes are recorded in the receipt. No system package was installed.

Activate the environment and build the required targets with explicit limits:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cmake --build /home/jay/build/nektar-v2026-mpi \
  --target IncNavierStokesSolver SolverUtilsUnitTests HITBoundaryUnitTests \
           HITBoundaryProbe HITFieldCompare --parallel 32
cmake --install /home/jay/build/nektar-v2026-mpi
ctest --test-dir /home/jay/build/nektar-v2026-mpi \
  -R '^(HITCoreUnitTests|HITBoundaryUnitTests|HITRestartUnitTests)$' \
  --parallel 1 --output-on-failure
```

The retained execution helpers run `base`, `refine` and `restart` phases,
configuration/source-checkpoint guards, and native norm extraction. Their
paths/hashes, actual argument vectors, fixture hashes, and named regression
commands are in the receipt. The committed generator/checker and native tools
provide the reusable test definitions. `--npz 1`, `--npz 4`, and `--npz 2`
exercise mesh, Fourier, and mixed layouts respectively at four ranks; eight-rank
mixed runs use `--npz 2`.

Compiler parallelism is 32; CTest concurrency is 1; solver runs use 1, 4 or 8
ranks with one OpenMP/BLAS thread. Overlapping validation drivers are bounded
by 18 scheduled solver workers in total. Builds/installations do not overlap
solver execution. These allocations remain below the hard 72-worker limit.

## Remaining scope

Milestones 5-8 remain open. This milestone establishes representation,
communicator, checkpoint and reproducibility behavior on straight-sided box
fixtures. Analytical Fourier solenoidality remains distinct from finite
polynomial projection and evolved solver divergence. Production spatial/time
resolution, domain/record-length effects, turbulence statistics and the
Comte-Bellot--Corrsin downstream-decay comparison remain later tasks. HDF5
checkpoint behavior and arbitrary curved-mesh serialization are not validated
by these runs.
