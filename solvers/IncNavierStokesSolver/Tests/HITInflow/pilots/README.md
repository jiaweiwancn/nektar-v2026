# Limited-band CBC inflow pilot tools

`prepare.py` writes one complete native affine quadrilateral session, an upstream
CBC1971 theta=42 spectrum copy and `manifest.json`. It does not run Nektar++.
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
wavenumbers rad/m. `--h` is the maximum affine element side length; actual
`hx`, `hy` follow integer element counts and are recorded. `--modes` is the
spectral/hp expansion size, so polynomial degree is `modes-1`. `--hom-modes`
is the even full real Fourier packing size, including DC and the reserved
Nyquist slot. Both transverse periods equal `--width`; y boundaries are native
reciprocal periodic conditions, and z uses canonical `HOMOGENEOUS=1D`.

The synthesis box has periods `(synthesis-length,width,width)`. The default
`--synthesis-counts auto` chooses even sizes with cutoff strictly below each
Nyquist wavenumber; an explicit `NX,NY,NZ` must meet the same constraint. The
solver z grid must also meet strict Nyquist. `--taper` is an absolute taper-start
wavenumber; default is `0.8*cutoff`. `REPEAT=Error` requires
`duration+2*dt<synthesis-length/12.7`. The declared synthetic record reserves
two timesteps beyond the job end to cover accumulated clock roundoff; sampling
uses the actual absolute solver time. Counts, seed, periods, spectrum, cutoff and taper
must stay fixed during solver resolution/time-step refinement. The continuum
mode count in the manifest is a cost estimate, not a discrete retained-mode
count. The inlet groups the retained Fourier sum by transverse wavenumber and
caches spatial phases; the reference direct sum remains available for checks.

The only experimental input is the committed upstream shell spectrum. Its
interpolation/extensions are `LogLinear`, `K4`, `Gaussian`, with `Preserve`
normalization. Mean velocity is 12.7 m/s and molecular viscosity 1.5e-5 m2/s.
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
Matrix` omits `USEFFT`, selecting the native matrix path.

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
