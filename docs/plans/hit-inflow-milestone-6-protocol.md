# Milestone 6 pilot protocol

Status: execution complete. The predeclared gates below are retained; the
[closure report](hit-inflow-milestone-6.md) records the failed 2% target and
explicitly limited production choice. This protocol is not the closure receipt.

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
   64. Keep generated fields and logs outside the source checkout.
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

## SVV pilots

Use independent spectral/hp and Fourier switches. The candidate uses the native
DG kernel in x-y with SVVDiffCoeff=1 and constant SVVVelocityMagnitude=(12.7,0,0).
This gives base viscosity 12.7*h/p times the native optimized modal kernel.
The z kernel is exponential, with separately recorded cutoff and viscosity.
Compare half/double SVV strength on identical realizations and sampling windows.
SVV sensitivity is a model uncertainty, not a physical-viscosity refinement.
The quadrilateral DG kernel uses a tabulated filter and ignores the numerical
SVVCutoffRatio value. Refining h/p or Fourier resolution also changes the SVV
operator, so resolution comparisons measure their combined effect. No claim of
a common fixed-viscosity PDE is attached to those comparisons.

The code review identified unequal damping of full-Fourier cosine/sine partners
when the kernel uses packed slot IDs. The candidate corrects that mapping to
physical wavenumber IDs and adds a sine/cosine analytic diffusion regression.
Require that regression before the SVV pilots. The flow diagnostic measures
native instantaneous SVV work separately in x-y and z, with a Fourier Parseval
check. This work is distinguished from temporal fluctuation dissipation.

All M6 jobs now share an affinity set of at most 64 logical CPUs, following the
user's stricter limit. The study runner permits at most 60 MPI ranks in total,
leaving capacity for its controller and serial postprocessing, with one OpenMP
and OpenBLAS thread per rank.

The benchmark also integrates the exact finite-lattice tensor multiplied by
1-sinc(U*kx*T/2)^2. This gives the expected continuous-time central variance for
a window T and explicitly accounts for time-independent kx=0 modes. Its bias
relative to ensemble variance is a finite-window effect, separate from the
conditional block sampling intervals. Avoiding exact recurrence alone does not
establish low-frequency convergence.

The kernel choice follows the local user guide and constant-Peclet design in
[Moura, Sherwin and Peiro (2016)](https://doi.org/10.1016/j.jcp.2015.12.009)
and the spatially developing flow study of
[Moura et al. (2020)](https://doi.org/10.1016/j.jcp.2019.109112).
These provide stabilization guidance, not a guarantee of CBC spectral fidelity.
The broader LES use of SVV is described by
[Karamanos and Karniadakis (2000)](https://doi.org/10.1006/jcph.2000.6552).

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
