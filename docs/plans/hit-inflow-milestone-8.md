# HIT inflow Milestone 8: documentation and reproducible examples

Status: **CLOSED — DOCUMENTATION AND REPRODUCIBILITY**, 2026-10-31.
The full 3D and Homogeneous1D examples, version-1 configuration reference,
source-data preparation, small-regression commands and experimental-comparison
reproduction workflow are delivered and exercised on Linux.

All eight implementation milestones are closed under their recorded scopes.
The [M7 comparison](hit-inflow-milestone-7.md) remains a descriptive limited-band
SVV result: full-band physical validation and 2% precision are unestablished.
The conditional narrow-band validation extension is deferred because the
full-band results do not justify it.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill / mode | academic-research-suite, experiment-agent / run and validate; PDF skill for rendered guide QA |
| Origin date / version | 2026-10-31 / hit_inflow_m8_v1 |
| Verification status | VERIFIED for the delivered examples, configured regressions and deterministic cached-data reproduction |
| Branch | `codex/hit-inflow-generator` |
| Tested code | `2e25b26b427bb63b77e630269798a34c770a4594` |
| Full-3D fixture source | `6f6c0f2e32666f5636aecdb9aa2c9c9f0a089f5a`; its XML bytes are unchanged in the final delivery |
| Final guide source | `4825defc36259ff2016d840a8995dba11bbd31b2` |
| Native solver/tools build | `29f467b6428ce5efc58d97d4cb481fbda0e432b1`; C++ unchanged, no rebuild/reinstall required |

The [machine-readable receipt](hit-inflow-milestone-8-validation.json) contains
exact commands, revisions, MPI/thread allocations, file hashes, outcomes and
output locations. Verification concerns these concrete checks. It does not
extend M7's physical-validation verdict or assert independent third-party CFD
reproduction.

## Entry points

| Deliverable | Location |
| --- | --- |
| Complete configuration reference | [User-guide HIT subsection](../user-guide/solvers/incompressible-ns.tex), label `sec:hit-inflow` |
| Full 3D von Karman example | [HITInflow3D](../../solvers/IncNavierStokesSolver/Examples/HITInflow3D/README.md) |
| Homogeneous1D tabulated example | [HITInflowQuasi3D](../../solvers/IncNavierStokesSolver/Examples/HITInflowQuasi3D/README.md) |
| Source data and small-regression instructions | [HIT test/data guide](../../solvers/IncNavierStokesSolver/Tests/HITInflow/README.md) |
| Pilot, spectrum and reproduction instructions | [Pilot guide](../../solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/README.md) |
| Four-case comparison driver | [reproduce_m7.py](../../solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/reproduce_m7.py) |
| Physical comparison and limits | [M7 report](hit-inflow-milestone-7.md), [M6 error budget](hit-inflow-milestone-6.md) |

The guide documents required attributes and units, Batchelor/von Karman/shell
spectrum conventions, normalization, tail extensions, strict Nyquist/time/record
bounds, solver restrictions, pressure/periodic boundaries, MPI transforms,
restart identity and the separate numerical/physical limits. The root and solver
example READMEs now point to the feature; older test READMEs use the current
64-worker execution limit.

## Small runnable examples

Both standalone sessions embed their meshes. The full-3D example has four
hexahedra; the Quasi-3D example has four xy quadrilaterals and eight Fourier slots
in z. They use degree 7, Ubar=1 m/s, nu=0.01 m2/s, dt=0.001 s and 20 steps,
mean-only initial volumes, and 12 physical history probes. The inlet sets
u/v/w with `HITInflow` and pressure with `H`; transverse boundaries are periodic,
and the outlet uses `HOutflow`. The synthetic box is `(2*pi,2*pi,2*pi)` m with
counts `(4,4,4)`, seed 20261028 and cutoff/taper 1.55/1.5 rad/m.
`TargetRMS=0.05 m/s` fixes expected total fluctuating energy. SVV is off.
These deliberately small examples are configuration demonstrations, not grid
experiment or turbulence-convergence cases.

The supplied full-3D session uses von Karman input. The supplied Quasi-3D table
is explicitly manufactured, with linear interpolation and zero extensions; it
contains no experimental values. Spectral/hp dealiasing is enabled in both.
Homogeneous dealiasing is disabled in the default Matrix session and the tested
FFTW comparison, so those transforms use the same numerical choices.

Five configurations completed t=0.02 s, each with 21 complete equal-time frames
(including t=0), finite u/v/w/p histories, finite native field diagnostics, and
verified evolved/final HIT metadata:

| Configuration | MPI ranks | npz | Backend | Final div RMS (1/s) | Relative flux imbalance |
| --- | --- | --- | --- | --- | --- |
| full3d-serial | 1 | 1 | 3D | 0.000395268 | 9.10558e-08 |
| full3d-mpi2 | 2 | 1 | 3D | 0.000395268 | 9.10559e-08 |
| h1d-serial | 1 | 1 | Matrix | 0.000427144 | 1.0499e-07 |
| h1d-mpi2 | 2 | 2 | Matrix | 0.000427144 | 1.0499e-07 |
| h1d-fftw-mpi2 | 2 | 2 | FFTW | 0.000427144 | 1.0499e-07 |

A serial `HITFlowDiagnostics` call reads each complete final field, including
parallel field directories. Maximum absolute u/v/w probe-history differences
were 1.3781e-11 m/s for full-3D serial/MPI, 2.2204e-16 m/s for H1D serial/MPI,
and 1.1416e-12 m/s for H1D Matrix/FFTW. All meet the preselected 1e-9 m/s small
consistency tolerance. These comparisons do not establish a physical resolution
or statistical uncertainty bound.

Two existing behavior details are now explicit:

- Native homogeneous `DEALIASING=True` requires `USEFFT=FFTW`. Matrix supports
  spectral/hp dealiasing alone. `prepare.py --backend Matrix --dealias` now fails
  before creating an invalid case, with a focused regression test.
- The legacy initial `_0.chk` at t=0 can precede HIT metadata initialization.
  An initial checkpoint without that metadata is not a verified HIT restart.
  The examples' evolved checkpoints at 0.01/0.02 s and final fields contain the
  correct identity and native build revision; retain all four fields and the
  stored absolute time when restarting.

The initial checker attempt and the rejected Matrix/dealiasing attempt remain
in `examples-01` and `examples-02`. Successful full-3D checks from the second
attempt were retained after verifying unchanged fixture/history hashes; the
three corrected H1D configurations ran in fresh `examples-03` folders. The
consolidated receipt records this provenance.

## Data, regression and analysis reproduction

On Linux, the documented workflow regenerated the seven CBC numeric fixtures
and manifest byte for byte, passed their read-only check and interface contract,
and reproduced the upstream-only von Karman fit exactly for its scientific
parameters. Table-2 longitudinal spectra remain reference data; Table-3
upstream shell data are the experimental inlet input. The exact Table-5
transcription is committed. No PDF download or new digitization is required.

All **13 Python tests** passed: eight reproduction-driver tests and five input
preparation tests. These exercise fresh-path and symlink protection, archive
immutability, explicit single-CPU execution, failure receipts/no retries,
fixed scientific arguments, and the invalid Matrix/dealiasing combination.

The configured CTest commands were exported into separate fresh working
directories and run sequentially. All **18 native CTest entries** passed:
14 solver regressions and four unit groups containing 48 cases. They cover
small seeded/tabulated, full-3D/H1D, MPI, transform, restart, manufactured-flow,
conservation, initial-field and SVV behavior. The exact commands and individual
logs are retained in `checks-01`; the test/data guide provides the same bounded
command-export workflow. No additional legacy suite or full production CFD
rerun was needed for these documentation/orchestration changes.

The new driver executed nine helper calls sequentially on CPU 50: four
analyses, four joins against the retained native profiles, and one five-figure
render. It reuses immutable NPZ caches by default and preserves all archive
inputs. New results live outside the original archive and source tree.

For all four cases, **14 scientific analysis subtrees and 15 native-work
subtrees match the original decoded JSON exactly**. All five figures' plotted
data also match exactly. New figure hashes match their new provenance. Paths,
commands and timestamps differ by design; complete SVG byte streams are not a
scientific equality criterion. The M7 qualification remains unchanged, including
32/32 failed sampling screens. No new fitting, independent ensemble inference
or improved physical agreement is claimed.

Run the documented complete reproduction command from any working directory:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
hit_src=/media/jay/新加卷/nektar-v2026
hit_archive=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-7/97014ab
hit_output=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/m7-reproduction-$(date -u +%Y%m%dT%H%M%SZ)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 "$hit_src/solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/reproduce_m7.py" \
  --run-root "$hit_archive" --output "$hit_output" --cpu 50 --source-root "$hit_src"
```

Reserve an available CPU in the inherited affinity. Output must be fresh and
outside the entire source/archive trees. The result contains per-case
`analysis.json` and `work.json`, five PNG/SVG pairs and plotting provenance,
logs and an atomic completed/failed receipt. The original archive uses absolute
paths; arbitrary relocation is not supported. Default-cache reproduction is
fully exercised. The optional `--rebuild-cache` orchestration path is unit-tested;
M8 did not repeat the 27.6 GB ASCII replay. The unchanged underlying cache reader
already created and validated every M7 production cache. Native late profiles
are reused in either mode, without new field solves or profiling.

The tested analysis environment is Python 3.10.12, NumPy 1.21.5, Matplotlib 3.5.1
and SciPy 1.8.0. Data conversion and driver orchestration use the standard
library; continuum fitting needs SciPy, analysis NumPy, and figures Matplotlib.
Native regression tools are built with `NEKTAR_BUILD_TESTS=ON`; the additional
unit groups require `NEKTAR_BUILD_UNIT_TESTS=ON`. The tested Release installation
has MPI and FFTW enabled. Source code is synchronized before checks; the native
C++ build remains the version verified in M6.

## User-guide build and rendered QA

The full user guide compiles successfully with the existing Linux LaTeX tools.
The final 337-page PDF and logs are in `user-guide-full-02`. Edited HIT pages
248--253 (printed pages 234--239) were rendered with Poppler and visually
inspected. The required spectrum-file header now fits on one literal monospace
line. Tables, XML listings, equations, page transitions and restart/example
instructions are readable, with no clipping or overlap in these pages.
The validation receipt records source/PDF/log/render hashes. The complete PDF
stays in Linux run storage, following the repository's generated-artifact policy.

The tested build command used one CPU:

```bash
hit_src=/media/jay/新加卷/nektar-v2026
hit_pdf=/media/jay/新加卷/nektar-v2026-runs/hit-inflow/user-guide-reproduction-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$hit_pdf/tmp"
cd "$hit_src/docs/user-guide"
TEXINPUTS="$hit_src/docs/user-guide//:" BIBINPUTS="$hit_src/docs:" TMPDIR="$hit_pdf/tmp" \
  taskset -c 52 latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir="$hit_pdf" user-guide.tex
```

## Execution and retained evidence

Authoritative host: `jay-Super-Server`, accessed using
`ssh jay@100.105.247.110`. The MPI Release environment remains
`/home/jay/.local/opt/nektar-v2026-mpi/env.sh`.

```text
M8 root: /media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-8/6f6c0f2
  checks-01/receipt.json                  Python/data/fit and native test commands
  checks-01/ctest-commands.json           Configured CTest entries
  examples-03/receipt.json                Five successful example configurations
  reproduced-m7-01/receipt.json           Sequential analysis/figure commands
  reproduction-comparison.json           Exact scientific-payload comparison
  user-guide-full-02/receipt.json         Final PDF/source/log/render bindings
  guide-page-QA-final/                    Edited-page PNGs
M7 input: /media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-7/97014ab
```

All raw fields, histories, NPZ caches, temporary files, logs and the compiled
manual remain outside the source checkout. Native regressions used at most eight
MPI ranks; examples at most two; comparison reproduction and guide compilation
one CPU each. Parallel numerical libraries used one thread. The maximum planned
simultaneous compute allocation across the validation stages was 13 workers,
with bounded controllers sharing headroom below the server's hard 64-worker
limit. No native installation was changed.

## Closure disposition

| M8 requirement | Evidence |
| --- | --- |
| Small runnable example per supported expansion | Full 3D and H1D delivered; all five serial/MPI/backend configurations pass |
| Configuration and boundary/restart reference | Expanded user-guide section compiles and passes rendered QA |
| Reproducible source preparation | Numeric fixtures/manifest byte-exact; upstream-only fit matches |
| Documented small regressions | 13 Python tests and 18 native CTest entries pass with exported commands |
| Regenerate experimental comparison | All four analyses/work payloads and five figures' data match exactly |
| Tested revision/settings and large-data locations | Recorded here and in the machine-readable receipt |
| Conditional narrow-band extension | Deferred because full-band criteria are unmet |

M8 completes the requested implementation/delivery sequence. Quantitative
physical validation still requires the upstream-only developed reference case,
better sampling qualification and the numerical controls recorded in M6/M7.
Those scientific requirements remain open.
