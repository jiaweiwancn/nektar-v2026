# Qualification rules for the limited-band HIT follow-up

Predeclared on 2026-10-29 during Milestone 6, before Milestone 7 experimental
scoring. Milestone 6 remains in progress until its long pilots, final build,
production manifest and error budget are complete. These rules supplement the
[pilot protocol](hit-inflow-milestone-6-protocol.md); they do not declare the
current candidate converged.

## Scope and numerical claims

The feasible candidate is an explicitly filtered, SVV-assisted calculation.
Preserve the upstream shell-spectrum amplitude and state the cutoff/taper,
retained lattice K/epsilon, physical viscosity and separate x-y/z SVV settings.
Missing inlet energy and dissipation are model limitations, not sampling error.
No full-DNS or full-spectrum CBC accuracy claim follows from this case.

Keep the provisional 2% target for a matched numerical refinement: each component
variance and K at every shared station, with the conditional intervals also
reported. Overlapping intervals alone do not demonstrate 2% equivalence.
State where a control was actually run. Near-field timestep/Fourier/buffer
sensitivities and a scalar IMEX model cannot bound the corresponding far-field
errors. Untested far-field effects must be labelled **not quantified**, rather
than assigned a zero contribution or combined into a numerical error bar.
DG h/p and homogeneous resolution changes also change the SVV operator.

Freeze the production input and these analysis rules before examining downstream
experimental discrepancies. Any later numerical or sampling refinement must be
motivated by its own diagnostics, with new inputs/provenance, rather than by
tuning to stations 98 or 171.

## Recording and sampling qualification

The candidate follow-up uses two independent seeds, 20261028 and 20261029, a
0.8 s startup and a subsequent 0.8 s recording interval. This duration is a
feasible selection to be tested at every station. The completed short-domain
records do not establish far-field stationarity, independent blocks, or a 2%
statistical uncertainty bound.

Before uncertainty-weighted experimental scoring:

1. Verify successful solver receipts, finite fields/histories, input hashes,
   complete equal-time records and absence of native iteration-limit warnings.
   Retain the declared recurrence margin, including startup and clock padding.
2. Compute each probe's own temporal-central variance and then average over the
   transverse plane. Never subtract an instantaneous plane mean. Report u/v/w
   variances and K, each seed separately, plus the between-seed range.
3. At every station, report full/half/last-75-percent window estimates and
   contiguous block results for 4, 8 and 16 blocks. Retain the full-window mean
   for block energies; also report own-block-centered estimates, excluded tail
   samples, block sizes, and observed block/temporal correlations. Distinguish
   the full-window variance from the mean block-energy estimator.
4. Flag a record as **not qualified at 2% sampling precision** if a component
   variance or K changes by more than 2% between the full and final-half window,
   or if a reported conditional 95% interval half-width exceeds 2% of the
   full-window statistic. These are engineering screening flags, not hypothesis
   tests. Passing them does not establish independent blocks or asymptotic
   confidence coverage. Report sensitivity to block choice even when a flag
   passes; two seeds do not justify a precise ensemble confidence interval.
5. If these diagnostics remain unstable, report discrepancies descriptively and
   retain conditional interval labels. Extend the record or independent
   ensemble, with a new recurrence-safe input when needed, before claiming
   precision unsupported by the available data. Do not silently remove a seed
   or choose a favorable time window.

Finite-window centering bias from the exact inlet tensor is a separate entry.
Its frozen kx=0 contribution and continuous-time expectation do not quantify
downstream sampling uncertainty. Spatial probe count is not an independent
sample count. Experimental uncertainty remains distinct; no universal
experimental error bar is prescribed here.

## Adjustment, spectra and dissipation

Use late native profiles and time-window diagnostics to state whether the
adjustment region is bounded, continues through the measurement domain, or is
indeterminate. Nominal convection transit alone does not establish startup
convergence. If adjustment is persistent or indeterminate, retain the option of
a development section whose reference-plane calibration uses upstream data
only; that would be a separately documented case.

Keep these observables separate: temporal-central TKE; instantaneous energy
about prescribed Ubar; molecular strain dissipation `2*nu*S:S`; native x-y/z
SVV work; and the z-spectrum derivative moment `sum(kz^2*variance_per_mode)`.
The last is not total dissipation. Instantaneous work is not a closed decay
budget without open-boundary transport, pressure and time-integration terms.

Label the direction and normalization of every spectrum/correlation. A discrete
z component spectrum is neither the three-dimensional shell input nor the
experimental longitudinal spectrum. Show retained-band limitations alongside
full-band experimental curves; they prevent a full-band validation verdict.
Milestone 7's full-band physical comparison remains pending or explicitly
scope-qualified, regardless of implementation-test success.
