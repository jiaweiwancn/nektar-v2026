# Milestone 6 pilot protocol

Status: execution in progress; this document is not a closure receipt.

Use upstream CBC1971 station 42 only to select inlet amplitudes and scales.
Do not calibrate parameters to stations 98 or 171 during this milestone.
Lengths, velocities, viscosity and angular wavenumbers are SI. Preserve the
tabulated spectrum amplitude, with LogLinear interpolation and explicit K4
and Gaussian tail policies. A smooth cutoff defines a filtered experiment;
it does not restore missing energy or dissipation.

## Gates and comparisons

1. Build on Linux, verify the grouped plane evaluator against the reference
   sum, and run the focused HIT and boundary regressions. Preserve absolute
   time, coefficient identity, conjugate pairing and restart metadata.
2. Measure generator construction, plane sampling, initialization and solver
   cost before allocating long runs. Check actual peak memory as well as mode
   storage. Record all MPI ranks and thread counts; total workers never exceed
   72. Keep generated fields and logs outside the source checkout.
3. Start with a short physical domain and cutoff near 100 rad/m. Independently
   vary inlet projection, spectral/hp resolution, homogeneous planes, timestep,
   transverse span and outlet buffer. Keep synthesis settings and seed fixed
   for discretization comparisons. A changed span or cutoff changes the
   realization and must be treated separately from discretization error.
4. Measure startup and spatial adjustment using initially consistent synthetic
   volume fields, with a mean-only initialization as a startup control where
   affordable. Use multiple time windows and nonoverlapping blocks. Establish
   the recording window empirically rather than assuming one eddy turnover
   time is sufficient.
5. Exercise a domain reaching the final comparison station, 8.322564 m from
   the reference plane, plus an explicit outlet buffer. Measure its resource
   use. Production duration must fit strictly inside the synthetic recurrence
   period, Lx/12.7, including startup; retain this inequality in every manifest.
6. Freeze the feasible production configuration and error budget before
   Milestone 7 comparisons. The provisional discretization target is less than
   2% change in each component variance and K, with overlapping sampling
   intervals. A failed target remains a measured limitation and requires either
   further refinement or an explicitly limited production claim.

## Statistics and uncertainty

At each sampling point subtract its own temporal mean before computing the
variance; then average the point variances over the transverse plane. Do not
subtract an instantaneous plane mean: that would erase valid zero-transverse
wavenumber fluctuations. Report component variances and K separately.
Use equal physical sampling times for paired discretization comparisons.
Report nonoverlapping block estimates and confidence intervals with the block
length, count, and observed sensitivity to doubling the block/window length.
Few blocks cannot establish an asymptotic uncertainty claim.

Examine one-dimensional spectral shape and resolved viscous dissipation
separately from K. Label spectra by their actual direction and transform;
do not identify a one-dimensional spectrum with the input three-dimensional
shell spectrum. Use native derivatives for volume dissipation/divergence.
Report raw divergence norms and a fixed physical normalization when comparing
different domains; volume-dependent normalization alone is insufficient.

Separate these error-budget entries: inlet projection; spatial, Fourier and
time discretization; transverse and outlet effects; startup/adjustment;
finite-record sampling; spectrum interpolation/tails and excluded scales;
stabilization/model dependence. A sampling confidence interval is not an error
bar for an unrepresented dissipation range.

## Scope if DNS is infeasible

The upstream Kolmogorov length is approximately 0.2905 mm for nu=1.5e-5 and
epsilon=0.474. Even k*eta=1 implies order one billion spatial points over a
9.3 by 0.6 by 0.6 m domain at uniform Nyquist spacing. This arithmetic is a
preflight estimate, not a measured solver benchmark. The final report must
combine the estimate with measured pilot memory/runtime.

For a limited-band configuration, state the retained angular-wavenumber band,
finite-lattice energy and dissipation, physical viscosity, all stabilization
settings (including absence of SVV), and unresolved energy/dissipation. Do not
claim a full DNS decay validation or fit extra viscosity to downstream data.
