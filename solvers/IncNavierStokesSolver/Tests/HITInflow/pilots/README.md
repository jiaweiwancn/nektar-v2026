# Limited-band HIT inflow pilot and comparison tools

`prepare.py` writes one complete native affine quadrilateral session, an upstream
CBC1971 theta=42 spectrum copy (or explicit von Karman parameters) and
`manifest.json`. It does not run Nektar++.
Use the Linux environment and run storage specified in the repository's
`AGENTS.md` for simulations and performance measurements. Run the solver from the
generated case directory so the relative spectrum filename resolves correctly.
The example below illustrates the CLI, not an accepted numerical configuration:

```bash
python3 prepare.py --output /path/to/new/run-directory \
  --length 1 --width 0.3 --h 0.1 --modes 6 --hom-modes 16 \
  --dt 0.0002 --duration 0.16 --record-start 0.08 --sample-every 10 \
  --seed 20261028 --cutoff 80 --taper 64 --synthesis-length 32 \
  --stations 0,0.5,0.8 --buffer 0.2 --probe-grid 8,16 \
  --backend FFTW --dealias --svv dg --svv-coeff 1 \
  --svv-cutoff 0.75 --svv-z-coeff 0.001 --svv-z-cutoff 0.75
```

All lengths are metres, times seconds, velocities m/s, and cutoff/taper
wavenumbers rad/m. `--h` is the default maximum affine element side length;
`--hx` and `--hy` independently override it on either axis. Actual
`hx`, `hy` follow integer element counts and are recorded. `--modes` is the
spectral/hp expansion size, so polynomial degree is `modes-1`. `--hom-modes`
is the even full real Fourier packing size, including DC and the reserved
Nyquist slot. Both transverse periods equal `--width`; y boundaries are native
reciprocal periodic conditions, and z uses canonical `HOMOGENEOUS=1D`.

The synthesis box has periods `(synthesis-length,width,width)`. The default
`--synthesis-counts auto` chooses even sizes with cutoff strictly below each
Nyquist wavenumber; an explicit `NX,NY,NZ` must meet the same constraint. The
solver z grid must also meet strict Nyquist. `--taper` is an absolute taper-start
wavenumber; default is `0.8*cutoff`. In addition to `REPEAT=Error`, preparation
requires the complete footprint `length+12.7*(duration+2*dt)<synthesis-length`.
The declared synthetic record reserves
two timesteps beyond the job end to cover accumulated clock roundoff; sampling
uses the actual absolute solver time. Counts, seed, periods, spectrum, cutoff and taper
must stay fixed during solver resolution/time-step refinement. The continuum
mode count in the manifest is a cost estimate, not a discrete retained-mode
count. The inlet groups the retained Fourier sum by transverse wavenumber and
caches spatial phases; the reference direct sum remains available for checks.
Wide planes can use bounded separable y/z phase tables within the same64MiB
auxiliary-cache cap, rather than allocating one index record per retained mode.

The CBC experimental inlet uses only the committed upstream shell spectrum.
Its interpolation/extensions are `LogLinear`, `K4`, `Gaussian`, with `Preserve`
normalization. The von Karman option also uses Preserve normalization with the
explicit upstream-only shape length and amplitude. Mean velocity is 12.7 m/s and molecular viscosity 1.5e-5 m2/s.
The finite cutoff/taper changes retained energy and dissipation. No downstream
data are loaded or used for tuning. A coarse limited-band run is not full DNS
of the CBC experiment and needs a documented numerical/model error budget.

Velocity correction uses nonlinear unsteady Navier--Stokes, IMEX2, Standard
extrapolation, continuous projection, and an iterative solve tolerance 1e-12.
`--pressure-solver XxtMultiLevelStaticCond` selects the tested direct pressure
option while retaining iterative velocity solves. `run_m6_case.py` then adds
the native `--use-scotch` reordering option. It limits its inherited CPU affinity
to 64 CPUs and fixes OpenMP/OpenBLAS to one thread per rank; concurrent jobs
must share the same total 64-CPU budget.
The inlet sets all three velocities with `HITInflow` and pressure with `H`;
the outlet uses native `HOutflow`. `--buffer` validates the distance from the
last requested station to the outlet and records it; it adds no sponge force.
`--dealias` enables both homogeneous and spectral/hp dealiasing. `--backend
Matrix` omits `USEFFT`, selecting the native matrix path. Native homogeneous
dealiasing requires FFTW, so `--backend Matrix --dealias` is rejected before
creating a case. Spectral/hp dealiasing alone is available in the solver for
matrix transforms, as demonstrated by the public Quasi-3D example.

SVV is disabled by default. `--svv dg` explicitly selects the spectral/hp
`DGKernel` and a separate homogeneous `ExpKernel`; no blanket SVV property is
used. The prescribed `SVVVelocityMagnitude=(12.7,0,0)` gives native constant
Peclet h/p scaling, multiplied by `--svv-coeff` (default 1). `--svv exp` selects
the spectral/hp `ExpKernel`: its coefficient is a direct viscosity, not that
DG scaling. Choose `--svv-coeff` explicitly for an exponential-kernel study;
the common CLI default 1 is not an experimentally calibrated viscosity.
The homogeneous exponential coefficient/cutoff are always explicit when SVV
is enabled (defaults 0.001 and 0.75). Every SVV knob is recorded even when off.
The quadrilateral DG kernel ignores the cutoff-ratio parameter: it uses its
native tabulated modal filter. Refining the resolution also changes this model.

`compare_study.py --root STUDY_DIR --reference CASE --output comparison.json`
compares completed case statistics, preserves their hashes and resource receipts,
and reports skipped failures. Paired block differences require the same canonical
inlet and aligned sampling windows; span, seed and window changes remain unpaired.
Its conditional intervals and 2% sensitivity flag are descriptive diagnostics.
They neither establish independent blocks nor close a validation milestone.

Initial conditions are the mean flow by default. `--initial-field FILE_OR_DIR`
writes a native `F VAR="u,v,w,p"` initial condition referencing an existing
field file or an MPI partition directory containing `Info.xml`. Paths and
individual file hashes are recorded; directory receipts sort relative filenames
and also hash their canonical JSON map (sorted keys, compact separators).
This schedule is designed for a
time-zero initial field. A restart file carrying nonzero `Time` needs a separate
review of absolute inflow/filter/sample times; this tool does not adjust those.

HistoryPoints writes `native.his` in physical velocity space at uniform
half-open yz probes, with x stations sorted and probe order station/y/z.
Output includes the initial state and every `sample-every` steps. Duration and
record-start must align with these steps. The default 8x8 probe grid need not
resolve every solver transverse mode; refine it independently for statistics.
Step checkpoints are disabled by default; `--checkpoint-every STEPS` enables
them for postprocessing. Time-based checkpoints remain disabled. The native
final solver field and cumulative Reynolds stress fields are retained.
Pass `--no-stresses` to omit cumulative full-field stress output explicitly;
the manifest records this choice and the final solver field remains enabled.

Native `ReynoldsStresses` uses `OutputStartTime=record-start`, nonmoving central
moments and `SampleFrequency=sample-every`. Its fixed output frequency is the
record-start step when positive. Thus the start snapshot contains one sample
and zero central variance; later snapshots accumulate samples since start.
The filter may also write intermediate snapshots at multiples of record-start,
and always writes its final field, possibly duplicating a periodic final
snapshot. The exact schedule is in the manifest and capped at 16 periodic
writes. With record-start zero only the native final snapshot is requested;
the first stress sample is at the first sampled solver step, while history also
contains t=0. Native stress and history moments should therefore only be compared
over identical sample sets. `OutputStartTime` gates sampling, not just output.

After a completed immutable solver run:

```bash
python3 analyze_history.py --manifest /path/to/run-directory/manifest.json \
  --output /path/to/run-directory/statistics.json --blocks 8 \
  --correlation-lags 1,2,4,8,16 --z-spectra
```

The analyzer streams the native history twice, validates variables, coordinates,
point ordering, finite values, frame completeness, uniform time steps and full
requested window coverage. A history path can be supplied as the positional
argument. `--start`/`--end` override the manifest's inclusive recording interval.
Statistics use each point's temporal mean and population variance (`M2/N`),
then average variances spatially; `K=0.5*(var_u+var_v+var_w)`. They never subtract
an instantaneous plane mean. JSON also contains spatial variation of point
statistics and first-half/second-half/last-75-percent startup diagnostics.

Blocks are equal-count, contiguous and nonoverlapping, with at least two
samples each. A trailing remainder is reported and excluded only from block
estimates. Block energies use the point's full-window temporal mean so that
block-specific mean removal does not suppress low-frequency energy; own-block
centered variances are also reported. Student-t 95% intervals describe the mean
of block energies **conditional on approximate stationarity and independent
blocks**. Nonoverlap does not establish independence. Inspect temporal and
block correlations, duration and block-length sensitivity before using them
in an error budget. They are not experimental uncertainty estimates.

Optional temporal correlations retain the fluctuating plane-mean mode and use
each probe's full-window temporal mean. Optional z spectra are one-sided
discrete variance per Fourier mode, include mode0 and a Parseval residual,
and require probe Nz at least the solver HomModesZ. They are z component
spectra, not spherical shell spectra or the experiment's longitudinal spectrum.
Their standard-library DFT costs O(Nz2) per sampled row. No automatic physical
acceptance thresholds are inferred. Provenance hashes cover prepared inputs,
the analyzer, history and manifest.

For dense completed records, `--numpy-transforms` explicitly selects optional
NumPy FFTs and vectorized temporal lag products. It retains the same streamed
validation, scalar moments, temporal centering, blocks and estimator definitions.
Its backend/version/dtype are recorded; floating-point reduction order may differ
for transforms and correlations. NumPy is imported only when requested. The
default analyzer continues to require only the Python standard library.

`plot_m6.py` renders measured short sensitivities, recording/block diagnostics
and completed long spatial comparisons as PNG/SVG with hashed provenance.
It requires Matplotlib. Use `--short-root`, `--sampling-root`, `--long-root`
and a fresh `--output` directory; `--skip-long` explicitly omits unfinished long
cases. It never turns a plotted sensitivity into a milestone decision.

`spectrum_budget.py --spectrum TABLE --output budget.json` uses SciPy quadrature
to compare continuum interpolation, low/high extensions, and the 100/200 rad/m
tapered bands without normalization. `--native-source-root SOURCE` adds hashes
of the native policy implementations when running an isolated script copy.
Its moments and derived isotropic length scales describe the continuum input
model; finite-lattice expectations, seeded records and solver results remain
separate quantities. Quadrature errors do not estimate input-model uncertainty.

`checkpoint_profiles.py CASE --checkpoints 2 4 6 --output NEW_DIRECTORY --cpu 7`
runs native `HITFlowDiagnostics` on completed partitioned checkpoints. It uses
one CPU and one OpenMP/OpenBLAS thread, copies session/spectrum inputs to
separate output directories, and hashes every checkpoint partition before and
after reading. A later checkpoint or a successful completed run receipt is
required so that the currently written checkpoint is not analyzed. It records
native build metadata, commands, flow JSON and complete x-column slab CSV
hashes. Slab energy is instantaneous about prescribed Ubar, rather than
temporal-central TKE. Use this with the active solver installation, before
replacing its libraries, to retain consistent diagnostic operators.

`plot_m6_diagnostics.py --long-root CASES --statistics-root SPECTRAL_JSONS
--profiles-root NATIVE_PROFILES --output NEW_DIRECTORY` adds a completed
long-domain PNG/SVG figure and hashed provenance. Each statistics file is named
`long-p9.json` or `long-p10.json`; each profile directory has the collector's
completed receipt for checkpoints at 0.85, 0.90 and 0.95 s. The tool checks input
identity/window, slab hashes and slab/global work closure. It plots temporal
TKE and discrete z mode energy separately from native molecular/SVV work and
fixed-reference divergence. Snapshot ranges are descriptive, not confidence
intervals. Run beside `plot_m6.py`, which provides shared validation helpers.

## Completed M7 models and reproducible comparison

The frozen production files are
[the CBC study](../../../../../docs/plans/hit-inflow-milestone-6-production-study.json)
and [the von Karman study](../../../../../docs/plans/hit-inflow-milestone-7-von-karman-study.json).
The [M7 protocol](../../../../../docs/plans/hit-inflow-milestone-7-protocol.md)
predeclares the windows, blocks, spectra and fixed-delay correlations. The
[report](../../../../../docs/plans/hit-inflow-milestone-7.md) gives all native,
MPI, thread and SVV settings and the qualified experimental comparison. Both
models decay too slowly; all 32 sampling screens fail the 2% threshold. M8
packages reproducibility of that result; it does not establish full-band
physical validation. Narrow-band validation remains deferred.

`prepare.py --spectrum von-karman --vk-length L --vk-amplitude A` selects the
analytic alternative. CBC remains the default. Both retain `Preserve`
normalization and the explicit energy cutoff/taper. To reproduce the upstream
fit (only theta=42 RMS and Lf) and numeric source preparation in fresh run
storage:

```bash
hit_src=/media/jay/新加卷/nektar-v2026
hit_test="$hit_src/solvers/IncNavierStokesSolver/Tests/HITInflow"
hit_work=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/source-preparation-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$hit_work"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 "$hit_test/prepare_cbc1971.py" --output-dir "$hit_work/cbc1971"
python3 "$hit_test/prepare_cbc1971.py" --check --output-dir "$hit_work/cbc1971"
python3 "$hit_test/check_contract.py"
python3 "$hit_test/pilots/fit_m7_upstream.py" --output "$hit_work/von-karman-fit.json"
```

Preparation/checks require only the Python standard library; the continuum fit
also uses SciPy. Analysis uses NumPy, and plotting uses Matplotlib with a
headless backend. The tested Linux versions are Python 3.10.12, NumPy 1.21.5,
Matplotlib 3.5.1 and SciPy 1.8.0. No PDF or downloaded data is needed: checked
source numbers, SI fixtures and the exact Table-5 transcription are committed.
The longitudinal Table-2 files remain downstream/reference data and cannot be
passed to the shell-spectrum inlet loader.

To regenerate **all four analyses, native-work joins and five PNG/SVG pairs**
from the existing authoritative archive, choose a fresh output directory outside
the source tree and outside the entire M7 archive:

```bash
hit_src=/media/jay/新加卷/nektar-v2026
hit_base=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-7/97014ab
hit_result=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/m7-reproduction-$(date -u +%Y%m%dT%H%M%SZ)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 "$hit_src/solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/reproduce_m7.py" \
  --run-root "$hit_base" --output "$hit_result" --cpu 50 --source-root "$hit_src"
```

CPU 50 must be available and in the inherited Linux affinity; change `--cpu`
if another allocation is required. The driver and its helpers use one CPU,
process the cases sequentially, and share the server's 64-worker limit with
other jobs. No solver or new native profile run is started. Default input is
four immutable `history.npz` caches and the already hashed three-snapshot native
profiles under `postprocessing-da607a9`. The original archive uses absolute
manifest/receipt paths; arbitrary relocation is not supported.

The output contains `<case-id>/analysis.json`, `<case-id>/work.json`,
`figures/` with five PNG/SVG pairs and `provenance.json`, `logs/`, and an atomic
`receipt.json` listing exact commands, environment, input/output hashes and
completed/failed steps. Failed attempts remain for inspection and are not
retried. Existing output paths are rejected before writing. Helpers verify
retained NPZ hashes; the driver checks large-file identities and bound small
inputs across the run.

Add `--rebuild-cache` to re-read the four native ASCII histories and build fresh
caches under the new output instead. This requires reading about 27.6 GB of
histories and storing about 5 GB of new caches, plus working memory. Keep that
output in Linux run storage. It still reuses the original late native profile
receipts; it does not regenerate the CFD fields. Without either the caches or
native histories and the native profiles, the four-run comparison cannot be
reconstructed solely from its compact committed report.

Reproduction means matching scientific arrays, statistics and plotted data.
New provenance has new timestamps, paths and command lines; SVG serialization
can also differ. Preserve each original result and its hash rather than
expecting a new output file's complete byte stream to match.


## CBC100 gated follow-up

The [100rad/m protocol](../../../../../docs/plans/hit-inflow-band100-protocol.md)
defines a short1m CBC pilot and the gates before any full-distance calculation.
It uses a600/480 radial input band, not a100rad/m radial cutoff. Its10m synthetic
period and fixed(2548,154,154) counts also support a matched800/640 comparison.
One case at a time shares the hard64-worker server budget with other jobs.

`run_band100.py` executes on Linux with an activated Nektar environment. It
requires a clean checkout on `codex/hit-inflow-generator`, exact full source
revision, explicit ranks/CPUs and a fresh output directory on the external disk.
For example, after native build/unit validation:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
hit_source=/media/jay/新加卷/nektar-v2026
hit_pilots=$hit_source/solvers/IncNavierStokesSolver/Tests/HITInflow/pilots
hit_revision=$(git -C "$hit_source" rev-parse HEAD)
hit_cpus=$(python3 -c 'print(",".join(map(str, range(48))))')
python3 "$hit_pilots/run_band100.py" --stage preflight \
  --run-root /media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/preflight-01 \
  --source-revision "$hit_revision" --ranks 48 --npz 16 --cpus "$hit_cpus"
```

The preflight runs a bounded2x2-point native generator benchmark, an inlet-projection diagnostic
and160 solver steps. It saves exact commands, executable/library hashes, per-rank
singleton CPU affinities and atomic running/failed/completed receipts. A failed
job is not retried. MPI ranks have one OpenMP/BLAS thread and share the selected
CPU pool with their bounded controller.
The controller audits all native main/helper threads. It records temporary
in-pool startup binding changes; after each `HIT_NATIVE_READY` readiness marker it
requires that rank's assigned singleton CPU. All ranks must be ready before
solver steps, and any CPU-pool escape fails immediately.
The global marker includes ranks without local inlet faces. `HIT_EVALUATOR`
reports only the inlet-owning subset and its cache methods.

After a completed preflight passes its12h forecast/180GiB sampled-memory gate,
run a separate fresh study using the same revision, binaries and allocation:

```bash
python3 "$hit_pilots/run_band100.py" --stage study \
  --run-root /media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/study-01 \
  --preflight /media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/preflight-01 \
  --source-revision "$hit_revision" --ranks 48 --npz 16 --cpus "$hit_cpus" \
  --cases base,dt
```

Every selected case receives its own fresh timing smoke before its0.6s recording.
Controls `x,z,bandwidth,svv,h` can be selected in order; `x` halves only streamwise
element spacing, and `bandwidth` needs the
common192-slot `z` reference. Separate studies can use `--reference-study` with
a matching completed study. A budget rejection stops the remaining expensive
recipe and preserves the evidence. Full-distance production is never automatic.

`--recipe original` preserves the original0.05m isotropic mesh/12.5us design.
Its completed native smoke forecast15.604h, so the12h gate held the recording.
The paired boundary-only test changed10--100rad/m band power by at most0.219%.
The separate `--recipe economical` keeps that inlet-plane resolution and tests
`hx=0.1m`, `hy=0.05m`, degree10 and25us. It requires a new preflight and matching
`--recipe economical` on the study command; original/economical preflights and
reference studies cannot be mixed. Its first study selects `--cases base,dt,x`.
Mesh and timestep savings are hypotheses until their native smoke and numerical
controls are inspected. The same12h/180GiB gates apply to every case.

Read progress without launching work:

```bash
python3 "$hit_pilots/run_band100.py" --status \
  /media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/study-01
```

The source helper `../check_band100_status.sh RUN_ROOT` runs this read-only
status command locally on Linux or through Tailscale SSH from the Mac.
`analyze_band100.py --manifest CASE/manifest.json --output FRESH_ANALYSIS_DIR`
uses sparse physical probes without a full transverse-grid requirement. It
reports all component temporal spectra, fixed10--100rad/m band powers, endpoint
values, sampling diagnostics and optional matched-control changes. Raw station42
Table2a is an upstream descriptive overlay only. Neither a short-pilot completion
nor a passed timing gate establishes downstream physical validation.

The pilot sets optional session parameters `HITReportEvaluator=1` and
`HITRequireAccelerated=1` to report cache methods and reject direct fallback.
`HITProbeDirectSamples=1` bounds direct references per boundary/rank/time in the
native diagnostic; projection errors compare raw/projected values at every
sampled boundary point. These execution/diagnostic controls do not alter the
realization identity or restore omitted energy.

For a boundary-only temporal projection check, use a full-duration prepared
session and set `HITProbeTemporalSamples=1601`, `HITProbeSampleStart=0.2`,
`HITProbeSampleDt=0.00025`, `HITProbeOutputPoints=16` and
`HITProbeDirectSamples=0`. Run `HITBoundaryProbe` with the same MPI allocation
and save its native receipt. Direct evaluation is skipped only after the
separate direct/cache preflight has passed. The probe checks finiteness and
projection maxima at all native points, exports a fixed dispersed subset, and
checks a repeated time separately from the uniform CSV grid. With a session
named `projection-temporal.xml`, analyze the fresh export with:

```bash
python3 "$hit_pilots/analyze_projection_band100.py" \
  --case /path/to/projection-case --stem projection-temporal \
  --output /path/to/fresh-projection-analysis
```

Inspect all four physical-band power changes before starting a CFD recording.
The 5% projection screen describes the exported probes; it does not establish
spatial convergence, experimental agreement or independent samples.

`review_band100.py` reviews the completed `base,dt,x` study without rerunning
CFD. It verifies command/input/log/cache hashes, reproduces the original moments
and spectra exactly, partitions temporal power, and adds explicitly separate
window and nested-probe sensitivity diagnostics plus PNG/SVG comparisons.
Use a fresh output directory outside the source and original study:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
python3 "$hit_pilots/review_band100.py" \
  --study /media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/f4b70ba/study-01 \
  --output /media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/f4b70ba/review-fresh \
  --source-root "$hit_source" --cpu 63
```

CPU63 must be available and share the64-worker server limit with other jobs.
The recorded helpers/protocol must retain their original hashes. The
[completed CBC100 review](../../../../../docs/plans/hit-inflow-band100-review.md)
holds full-distance production: the limited spectral screens pass, but the
energy/refinement and sampling gates are not closed. Nested8/16/32 probe grids
cannot establish convergence to an unrecorded64/128 grid.
