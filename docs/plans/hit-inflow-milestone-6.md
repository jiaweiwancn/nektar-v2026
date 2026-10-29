# HIT inflow Milestone 6: resolution, domain and resource study

**Status: IN PROGRESS.** The full-distance pilots must finish before the
production configuration and Milestone 6 closure can be signed off. This report
records completed evidence without treating timing trials or analytical models
as downstream solver validation.

This implements Milestone 6 of the [implementation plan](hit-inflow-implementation-plan.md)
under the [predeclared protocol](hit-inflow-milestone-6-protocol.md). No downstream
experimental values have been fitted or scored. Milestone 7 remains separate.

## Scope and provenance

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Mode | implementation and numerical verification |
| Branch | `codex/hit-inflow-generator` |
| C++ build and 47 unit cases | `d74cae8` |
| Corrected 22 registered solver regressions | `812a2c8` |
| Long and follow-up solver inputs/runners | `812a2c8` |
| Optional spectral analysis | `a3caa84` |
| Authoritative host | `jay-Super-Server`, accessed through `ssh linux-usb` |
| Main run storage | `/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-6` |

The C++ sources are unchanged between `d74cae8` and `812a2c8`. Later changes in
`a3caa84` concern optional Python postprocessing. The active solver checkout
remains clean at `812a2c8`; the optional analyzer was copied into external run
storage with its source hash, without replacing running solver libraries.

All M6 work uses at most **64 CPUs**, following the user's stricter limit.
Builds use 32 jobs and never overlap solvers. OpenMP and OpenBLAS use one thread
per MPI rank. The current allocation is two long cases of 28 ranks each
(`npz=2`) plus one short case of four ranks (`npz=4`): 60 MPI ranks, with four
CPUs reserved for controllers or bounded postprocessing. A saved
`live-affinity-change.json` records placement across both sockets after the
initial common CPU set was found to occupy one socket. Each long rank is pinned
to one physical core; the combined allowed set still contains exactly 64 CPUs.
No session, field or solver binary was changed by that placement adjustment.

## Feasible physical scope

The upstream-only specification uses mean speed 12.7 m/s, molecular viscosity
`1.5e-5 m2/s`, the committed CBC1971 station-42 shell spectrum, log-log
interpolation (`LogLinear`), `K4`/Gaussian extensions, and `Preserve`
normalization. The synthesis periods are `(32,0.6,0.6) m`; the cutoff/taper
pair is `(100,80) rad/m` in the current long pilots.

| Quantity | Cutoff 100, taper 80 | Cutoff 200, taper 160 |
| --- | ---: | ---: |
| Retained independent modes | 97,807 | 777,106 |
| Native mode payload per generator | 8,607,016 bytes | 68,385,328 bytes |
| Expected periodic K, m2/s2 | 0.02599846997 | 0.04393868581 |
| Realized periodic K, m2/s2 | 0.02608822850 | 0.04394070614 |
| Expected molecular dissipation, m2/s3 | 0.00265132404 | 0.01198565572 |

These are native finite-lattice generator measurements. With the upstream
reference `K=0.073926 m2/s2` and `epsilon=0.474 m2/s3`, the 100-rad/m choice
retains only about **35.2% of energy and 0.56% of dissipation**; the 200-rad/m
choice retains about 59.4% and 2.53%. The reference totals and finite-lattice
integrals are different quantities; these percentages describe the severe
represented-band limitation, not a measured downstream error bar. No retained
band is silently rescaled to restore the missing energy.

The upstream Kolmogorov length is approximately 0.2905 mm. Even the planning
threshold `k*eta=1` would require roughly 3.96 billion retained modes in this
synthesis box, or about 349 GB of mode payload **per rank**, exceeding the
server's approximately 251 GiB RAM before fields and matrices. Uniform Nyquist
spacing over the 9 by 0.6 by 0.6 m domain gives approximately 4.4 billion points.
These are arithmetic planning estimates, not allocations or measured DNS runs.
They motivate a limited-band study with explicit SVV, rather than a full DNS
claim. SVV stabilizes represented modes; it does not supply the omitted inlet
energy or guarantee the experimental cascade.

## Implementation and verification

The inlet evaluator now groups modes by transverse wavenumber, retains distinct
streamwise temporal phases, and caches spatial phases within a 64 MiB payload
limit. Evaluation uses absolute time and preserves coefficients, random-number
identity, conjugate pairing and restart metadata. Unsupported/cache-limited
cases retain the direct-sum fallback. The native constrained boundary projection
is unchanged.

On Linux, a 512-point timing check measured 1.143 ms for the grouped plane and
6.425 ms per point for the original full reference evaluation. Extrapolating the
sparse direct check to 512 points gives approximately 2,880 times the grouped
cost. This is a plane-sampling comparison, not a solver speedup. Maximum velocity
differences were below `5.1e-13 m/s` and variance differences below `1.7e-15 m2/s2`.

The full-Fourier SVV kernel now assigns the same damping to each cosine/sine
pair, using its physical wavenumber rather than its packed slot ID. An analytical
phase-pair diffusion regression checks this behavior. Reduced Fourier bases keep
their existing convention. The native DG kernel and Fourier work diagnostics
are checked against the assembled stiffness quadratic form and Parseval's
relation. Native initialization also now correctly unpacks assembled global
coefficients before imposing local Dirichlet values.

All **47 unit cases** pass in four suites, and all **22 selected solver
regressions** pass, including eight legacy controls. The regression batch uses
CTest with eight worker slots and exact registered Tester commands in external
working directories. Two initial failures exposed checking issues: the first
absolute-norm fixture used a parameter override that left the parsed mean
expression unchanged; the constant-mean MPI pressure error was `8.34e-12`, above
the old signed-maximum bound. An explicit second-mean XML fixture fixes the
former; the latter now has a `2e-11` pressure bound, while velocity bounds remain
`5e-12`. Helpers now compute maximum absolute errors over the full communicator.
Earlier signed-maximum initialization diagnostics are excluded from Linf claims.

The explicit NumPy spectral backend agrees with the standard-library reference
on a native history and independent even/odd transform fixtures. Scalar means,
variances, blocks and windows are identical. Transform/correlation differences
are below `1.1e-15` in absolute value; maximum native Parseval residual is
`4.5e-16`. Both paths retain each probe's own temporal mean and preserve the
zero-transverse-wavenumber contribution.

## Completed short-domain sensitivities

The baseline has length/span 0.6 m, maximum element side 0.1 m, polynomial
degree 7, 32 full Fourier slots, `dt=1e-4 s`, and a 0.1--0.2 s record. It uses
the same 100-rad/m synthesis and seed for discretization controls. The following
are maximum changes in any component variance or K over x=0, 0.15 and 0.45 m.
Each isolated control is compared with its own baseline; comparing every row
with the finest mesh would confound the changed parameter with spatial error.

| Isolated change | Maximum relative change |
| --- | ---: |
| Degree 7 to degree 9 | 2.0267% |
| Degree-9 result versus h=0.05 m, degree-7 result | 0.1695% |
| Baseline timestep halved | 0.5269% |
| Fourier slots doubled | 0.00170% |
| Domain length increased to 0.9 / 1.2 m | 0.00296% / 0.00313% |
| Mean-only instead of Gaussian-consistent initialization | 0.0000177% |
| SVV coefficients halved / doubled | 0.2515% / 0.2537% |

The degree-7 baseline fails the provisional 2% spatial target: its w variance
at x=0.45 m is 2.1245% below the h-refined reference. Degree 9 is a better
candidate. Refining h/p or Fourier resolution also changes the SVV operator;
these are coupled numerical/model sensitivities.

Inlet projection was measured separately with native quadrature points at five
absolute-time calls, including repeated/reversed calls. The raw direct/grouped
agreement is below `6.7e-13 m/s`, with identical repeated coefficients and
nonzero Fourier content in all velocities. Maximum projected errors are:

| Representation | Maximum over velocity components, m/s |
| --- | ---: |
| Baseline degree 7 | 0.0010390 |
| Degree 9 | 0.00005707 |
| h=0.05 m, degree 7 | 0.000006562 |

Integrated-energy agreement does not establish dissipation agreement. At the
final short snapshot, total instantaneous SVV work divided by molecular strain
dissipation is approximately 4.679, 0.1601 and 0.03766 for the baseline,
degree-9 and h-refined cases. The corresponding fixed-reference divergence norms
are approximately 0.009395, 0.001868 and 0.001233. These distinguish converging
spatial representation from exact pointwise discrete solenoidality. Small net
flux/Gauss residuals alone do not establish the latter.

The no-SVV control failed with NaNs and repeated native iteration-limit warnings.
It is retained as a rejected stability control; it establishes no universal
claim that every unstabilized configuration is unstable. Span changes alter the
mode lattice and realization, even at the same seed. The first 0.9/1.2 m span
results differ from the fine reference by up to 8.25%/6.93%; they are unpaired
box-plus-realization sensitivities. Second-seed controls are still being run.

## Recording length and uncertainty

Statistics use each fixed probe's temporal-central variance, then a transverse
average. Instantaneous plane means are never removed. The extended short case
reproduces its initial 0.1--0.2 s statistics across the changed MPI layout and
allows longer windows without changing the inlet realization.

| Record at x=0.45 m | K, m2/s2 | Samples |
| --- | ---: | ---: |
| 0.1--0.2 s | 0.02409977658 | 201 |
| 0.1--0.4 s | 0.02502467584 | 601 |
| 0.1--0.8 s | 0.02545261927 | 1,401 |

The 0.1 s record is too short to represent the seeded ensemble. For the 0.7 s
record, eight-block conditional component interval half-widths at x=0.45 m are
2.24%, 2.08% and 2.77%. Four blocks give 2.54%, 2.35% and 4.53%; sixteen blocks
give 2.01%, 2.27% and 2.09%. Some block correlations remain substantial. These
intervals are conditional diagnostics, not established independent-sample
confidence bounds or 2% equivalence tests.

The exact finite-lattice tensor also predicts continuous-time centering bias:
for a 0.1 s inlet window the expected component biases are -7.68%, -3.94%,
-3.94%; for 0.8 s they are -0.983%, -0.493%, -0.493%. Frozen kx=0 modes alone
contribute approximately 0.313%/0.156%/0.156% of ensemble covariance. This
analytical expectation differs from both the actual seeded record and sampling
intervals. Avoiding recurrence does not by itself establish low-frequency
convergence.

## Full-distance work still required

The current 9 m pilots use degrees 9 and 10, h=0.1 m, 32 Fourier slots,
`dt=2.5e-5 s`, both dealiasing switches, and an absolute 0.8--0.95 s record.
They reach x=3.612896 and 8.322564 m with a 0.677436 m outlet buffer and retain
step checkpoints every 0.05 s. The 0.8 s start exceeds the 0.70866 s nominal
domain transit; actual startup and spatial adjustment still require inspection.

The smaller timestep follows a scalar IMEX2 calibration using the native
lattice/spectrum, not a replacement for Navier--Stokes validation. At the final
convection age, `dt=5e-5` gives modeled transverse variance error about 1.40%
but correlation about 0.9595; `dt=2.5e-5` reduces these to 0.179% and 0.9974.
Space, nonlinear evolution, pressure, SVV and boundary effects are excluded
from that analytical model. The fine short-case timestep chain and matched
long-degree comparison must provide the corresponding solver evidence.

Native checkpoint postprocessing will provide complete x-column slabs of
molecular strain dissipation, separate x-y/z SVV work, divergence and
instantaneous energy about prescribed Ubar. That energy is explicitly labelled
and cannot be substituted for temporal-central TKE. Direction-labelled z
spectra and correlations are computed separately from immutable histories.

The long pilots and matched linear-solver resource trial are still running.
The final production choice, measured runtime/RAM/storage budget, adjustment
assessment, uncertainty criteria and machine-readable closure receipt remain
outstanding. A proposed longer production record must stay inside the
`32/12.7 = 2.519685 s` synthesis period, including startup and clock margin.

## SVV interpretation

The candidate follows the local incompressible-solver guide's recommended
spectral/hp DG kernel, with prescribed speed `(12.7,0,0)` and coefficient 1.
The native base viscosity is `12.7*h/p` times the optimized modal filter. The
quadrilateral table, rather than the numerical cutoff-ratio setting, defines
that filter. The homogeneous exponential coefficient is independently
`0.0238125 m2/s` with ratio 0.75; it is not multiplied by the DG h/p scaling.

The constant-Peclet/optimized-kernel rationale is supported by
[Moura, Sherwin and Peiro (2016)](https://doi.org/10.1016/j.jcp.2015.12.009)
and [Moura et al. (2020)](https://doi.org/10.1016/j.jcp.2019.109112).
The latter's spatially developing model grid flow does not establish fidelity
for this three-dimensional CBC case. The broader SVV LES framework is described
by [Karamanos and Karniadakis (2000)](https://doi.org/10.1006/jcph.2000.6552).
Measured instantaneous SVV work is a model diagnostic; open-boundary transport,
pressure, time integration and molecular dissipation remain distinct terms in
the energy budget.
