# HIT inflow Milestone 7: downstream experimental comparison

Status: IN PROGRESS. The CBC production pair is running; the frozen von Karman
pair is queued after its successful completion. This file is not a closure
receipt. Follow the [predeclared protocol](hit-inflow-milestone-7-protocol.md)
and [M6 qualification rules](hit-inflow-milestone-6-qualification.md).

## Scope and experimental evidence

The calculation is a limited-band SVV study. M6's full-distance degree refinement
changed a component variance by 4.154%, exceeding the 2% target. Far-field
timestep, Fourier/SVV cutoff, span, buffer and SVV-strength errors remain
unquantified. Only about 35.2% of the upstream experimental K and 0.56% of its
dissipation enter the CBC case. Full-band experimental overlays therefore
describe discrepancies; implementation-test success does not establish physical
fidelity.

The principal experiment is Comte-Bellot and Corrsin (1971), the 5.08 cm grid
with upstream age theta=42. Table 2(a) contains the one-sided longitudinal
spectrum; Table 3 contains the isotropically reconstructed shell spectrum;
Table 4 contains RMS, decay-derived dissipation and inferred scales. The
[original transcription](../../solvers/IncNavierStokesSolver/Tests/HITInflow/cbc1971-source.json)
was rechecked against complete PDF pages. The separately saved
[Table-5 full-band correlation](../../solvers/IncNavierStokesSolver/Tests/HITInflow/cbc1971_full_band_correlation.json)
was transcribed from printed page 302, PDF page 30. No figure digitization or
universal experimental error bar is used.

The measurement-region speed is 12.7 m/s; approach speed U0=10 m/s defines
experimental age. Thus theta=42+U0*x/(M*Ubar), with M=0.0508 m. The downstream
comparison planes are x=3.612896/8.322564 m for theta=98/171. Table 5 instead
uses xi=x/M, giving xi=71.12/163.83 at these planes. Any interpolated experimental
correlation is labelled as interpolation rather than a published value.

## Frozen model comparison

Both families preserve amplitudes and apply the same radial energy cutoff/taper
of 100/80 rad/m. Both use the M6 degree-10, 32-slot Fourier, 9 by 0.6 by 0.6 m
configuration, molecular viscosity 1.5e-5 m2/s, IMEX2 dt=25 microseconds,
both dealiasing switches, DG coefficient 1 and independent homogeneous
exponential viscosity 0.0238125 m2/s with cutoff ratio 0.75.

The CBC inlet uses only the station-42 Table-3 shell spectrum, log-log
interpolation, K4/Gaussian extensions and Preserve normalization. The
[von Karman fit](hit-inflow-milestone-7-upstream-fit.json) uses only upstream
u_rms=0.222 m/s and Lf=0.024 m. Its native shape length is 0.03213564670 m and
amplitude is 0.7083410722 m^(4/3)/s2. The uncut model has K=0.073926 m2/s2 and
Lf=0.024 m; its inertial-tail dissipation diverges. The common energy taper
regularizes the simulated band without restoring omitted energy.

| Input diagnostic | CBC | von Karman |
| --- | ---: | ---: |
| Native retained modes | 97,807 | 97,807 |
| Native expected periodic K, m2/s2 | 0.02599846997 | 0.02494910477 |
| Native expected molecular dissipation, m2/s3 | 0.00265132404 | 0.00251070503 |
| Continuum filtered longitudinal integral scale, m | 0.050112 | 0.051813 |

These are prescribed input-model/lattice quantities, not evolved flow results.
The native VK parameter/evaluator preflight passed with maximum cached/direct
velocity difference 2.683e-13 m/s. The four new preparation contract tests pass
on Linux. No C++ source or installed binary was changed for these Python input
and analysis extensions.

Use the two disjoint input-counter streams 20261028 and 4315228324 in each
family. This avoids direct counter reuse between streams without proving
statistical independence. Matching streams across the two model families are
common random inputs for the comparison; four cases are not four independent
replicas.

## Execution provenance

The CBC pair was launched at 2026-10-29 12:58:07 UTC from source `97014ab`, with
native binaries built at `29f467b`. The intervening commit is documentation only.
The upstream fit/protocol and VK preparation/recipe were committed as `933ace5`
and `ca71712`, before any new downstream experimental scoring.

Authoritative run storage is
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-7/97014ab`.
The host is accessed over Tailscale with `ssh jay@100.105.247.110`.
Each production pair has two 24-rank jobs, npz=8 and OMP/OpenBLAS one. Every rank
and its helper threads are pinned to one physical core. The common 64-CPU set
covers both sockets; bounded analysis and other observed server jobs share the
remaining worker budget. Exact launch commands, wrapper hashes, native
metadata, affinities and the concurrent-job audit are retained in run storage.

Both jobs run from time zero to 1.6 s, with primary statistics over 0.8--1.6 s.
The declared synthetic interval and full space-time footprint remain below the
32 m repetition limits. Histories contain 64 by 32 transverse probes at each
of eight stations, sampled every 0.0005 s. Raw data, fields, checkpoints and
large analysis caches remain outside the source checkout.

## Analysis implementation checks

At source `da607a9`, all 15 focused analyzer/native-work tests pass on Linux
(eight analyzer and seven collector tests), alongside four preparation tests.
They check one-sided even/odd FFT normalization and the frequency/wavenumber
Jacobian, periodic-Hann energy normalization, exact frozen-advection lag sign,
window/block flags, immutable-input protection, periodic correlation conventions,
analytical native slab closure, interpolation and derived-scale units.

The completed M6 degree-10 native history was independently processed through
the new cache/analyzer. Its 301 selected samples reproduce the saved component
variances and K within 6.94e-18 m2/s2; the maximum boxcar Parseval residual is
2.22e-16. This pipeline check uses an explicit 128-sample Welch segment for the
short fixture; production retains the frozen 512-sample choice.

A matched full-field comparison at t=0.1 s gives exactly zero L2 and maximum
absolute differences in u/v/w/p between the original `d74cae8` pilot and the
new `29f467b` production build. Checkpoint hashes were unchanged during reading.
The declared recording durations differ, so this is a field/control comparison,
not a claim that the restart identity metadata is interchangeable.

A separately labelled synthetic analytic-signal fixture exercises all four
model/stream labels and renders the comparison plots. Its values are not CFD
results and are excluded from production scoring. Figure layout and source
conventions were checked before production data become available. Native work
plots and the actual production figures still require final inspection.

## Analysis and qualification

The production comparisons, sampling diagnostics and native late profiles
remain to be collected before closure. Required outputs are per-stream
full/final-half/last-75-percent statistics, 4/8/16-block sensitivity, conditional
intervals, PSD/Parseval checks, represented scales, absolute and normalized
decay, Table-2 spectral overlays and Table-5 delayed correlations. Molecular
strain dissipation and x-y/z SVV work remain separate from temporal-central K.

The far station's delay is 0.65532 s. Requiring both signals within 0.8--1.6 s
leaves only 0.14468 s of paired overlap. An additional explicitly labelled
forced-inlet extension may use earlier recorded inlet samples with post-start
downstream samples; this does not establish downstream stationarity. Finite
periodic correlation integrals, isotropy/Taylor assumptions and unclosed energy
transport terms must be disclosed separately.

A quantitative physical-validation claim requires an upstream-only developed
reference-plane alternative if adjustment persists. An unfavorable, unqualified
comparison can be closed as a reproducible descriptive result with its limits;
it cannot be converted into agreement by downstream fitting. Milestone 8 remains
separate.
