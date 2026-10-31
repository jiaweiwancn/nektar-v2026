# HIT inflow Milestone 7 comparison protocol

Status: COMPLETED, 2026-10-31. The [M7 report](hit-inflow-milestone-7.md) and
[validation receipt](hit-inflow-milestone-7-validation.json) close the prescribed
work as a descriptive limited-band comparison; full-band validation and 2%
precision are not established. The choices below remain the predeclared record,
written before examining the new production runs against downstream values.
It continues the [implementation plan](hit-inflow-implementation-plan.md) and
retains the [M6 qualification rules](hit-inflow-milestone-6-qualification.md).

## Inputs and execution

Use the frozen [CBC recipe](hit-inflow-milestone-6-production-study.json) without
retuning: 9 m domain, 0.6 m periodic spans, degree 10, 32 Fourier slots,
IMEX2 dt=25 microseconds, molecular viscosity 1.5e-5 m2/s, mean speed 12.7 m/s,
100/80 rad/m radial energy cutoff/taper, DG coefficient 1 and separate Fourier
exponential viscosity 0.0238125 m2/s with cutoff ratio 0.75. Preserve upstream
amplitudes. The source revision at CBC launch is `97014ab`, with native binaries
built at `29f467b`; the intervening commit contains documentation/assets only.

Run seeds 20261028 and 4315228324, retaining each separately. Their enclosing
input-counter ranges are disjoint; this avoids direct counter reuse without
proving statistical independence. Record the full 0--1.6 s history, and use
0.8--1.6 s for the primary station statistics. Sampling is every 0.0005 s on
64 by 32 half-open transverse grids at the eight prescribed stations.

The CBC pair starts first. Follow with a von Karman pair using the same mesh,
solver, viscosity, cutoff/taper, duration, probes and two streams. Fit its
uncut model only to station-42 u_rms=0.222 m/s and Lf=0.024 m. For
E(k)=A*L^(5/3)*(L*k)^4/[1+(L*k)^2]^(17/6), define
I4=Beta(5/2,1/3)/2 and I3=Beta(2,5/6)/2. Set
Lf=(3*pi/4)*L*I3/I4, K=3*u_rms^2/2 and A=K/(L^(2/3)*I4).
Apply the common radial taper afterwards, without restoring omitted energy.
The uncut inertial tail has divergent molecular dissipation; the represented
band's finite dissipation must be measured separately. No downstream values or
experimental decay law enter this fit.

Use Linux through `ssh jay@100.105.247.110`. Each production pair has two
24-rank jobs, `npz=8`, OMP/OpenBLAS one. Pin each rank and its helper threads to
one physical CPU, with 48 solver cores and a common 64-CPU allocation covering
controllers and bounded analysis. Store large output under
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-7/`.
Preserve failed/interrupted attempts, exact commands, input/output hashes and
native build metadata. Do not replace the installation during live CFD jobs.

## Quantities and sampling

Use each probe's own temporal mean, then average over each transverse plane.
Report u/v/w variances, RMS, K, mean velocities, Reynolds covariances and
anisotropy, separately for each stream. Retain full/final-half/last-75-percent
window estimates and 4/8/16-block estimates, own-block centering, unused tails,
conditional intervals and observed correlations. Apply the frozen 2% screens
at every station. Spatial probes and two pseudorandom streams are not a known
independent sample count. Failed precision screens require descriptive
comparisons and explicitly conditional intervals, not a fitted uncertainty bar.

Compare absolute Table-4 RMS and K, normalized decay, and inferred decay rates.
Keep measured molecular strain dissipation, x-y/z SVV work, temporal-central K
and energy about prescribed mean separate. Derived Taylor/Reynolds/Kolmogorov
scales must identify the dissipation estimator used. A temporal Taylor-hypothesis
estimate or z-only derivative moment is not a native full-gradient measurement.
Finite periodic correlation integrals must state the integration limit and
zero-mode treatment; do not silently equate them to the experiment's inferred
unbounded-field integral scales.

Use late native checkpoint profiles to assess adjustment and spatial work.
A nominal domain flush does not prove equilibrium. Report the existing 4.154%
long-domain refinement difference, coupled DG/model change, and unquantified
far-field timestep, Fourier/SVV-cutoff, span, buffer and SVV-strength effects.
Do not add these dependent effects into a total numerical error bar.

## Spectra and convection-delayed correlation

Estimate one-sided temporal PSDs of local-mean-subtracted u, then average over
the plane. Use a rectangular full-record periodogram for exact discrete
Parseval closure, and a Hann-window Welch estimate with 512 samples per
segment and 50% overlap for the primary spectral overlay. Record segment count,
frequency spacing, window-energy normalization and Nyquist. Convert to the
experimental one-sided longitudinal convention using k1=2*pi*f/U and
F11(k1)=U*Sf(f)/(2*pi). The saved histories have no explicit anti-aliasing filter
before every-20-step sampling; retain that limitation and the finite Nyquist
when plotting the full experimental Table-2 curves.

The homogeneous z spectrum is a different observable. A low-k1 restriction of
the published longitudinal spectrum is not a shell-filtered reference, because
F11(q)=integral_q^infinity E(k)*(1-q^2/k^2)/k dk. Any optional filtered reference
must taper the shell spectrum first and state isotropy/continuation assumptions.
Filtering downstream experimental data also differs physically from evolving
a filtered inlet with SVV. All full-band overlays remain descriptive.

Use the principal 5.08 cm grid's exact Table-5 transcription. Its abscissa is
xi=delta_x/M=U*delay/M, whereas experimental age is theta=U0*t/M with U0=10 m/s.
At identical transverse coordinates, pair the upstream signal at t with the
downstream signal at t+delta_x/U. Linearly interpolate fractional-sample delays,
report the convention, paired sample count, interval and effective overlap.
For the primary correlation, require both observations inside the declared
post-start record. At x=8.322564 m the delay is 0.65532 s, leaving only
0.14468 s of common post-start overlap. Never label this an 0.8 s estimate.

Also report a separately labelled extension pairing post-start downstream
signals with earlier recorded forced-inlet signals. Stationarity by construction
of the imposed inlet does not establish downstream stationarity or eliminate
finite-window bias. Preserve own-pair-window means and variances and distinguish
that Pearson estimate from full-window normalization. Report interpolation of
Table-5 experimental values explicitly; the station-98/171 positions are not
published Table-5 abscissae. No narrow-band validation is added in this milestone.

## Interpretation and closure

M6 established neither full-spectrum fidelity nor 2% numerical convergence.
The new runs test the frozen feasible configuration and its statistical
qualification. If the adjustment region persists, an upstream-only calibrated
development-region case is required before a quantitative physical-validation
claim. An unqualified or failed comparison is saved as a descriptive result;
it does not justify tuning to downstream stations or asserting validation.

Close M7 only after both prescribed spectrum families, reproducible Table-2,
Table-4 and full-band Table-5 comparisons, sampling diagnostics, represented
scales, separate molecular/SVV work and unresolved discrepancies are saved.
Keep any physical-validation verdict scope-qualified. M8 documentation/example
delivery remains separate.
