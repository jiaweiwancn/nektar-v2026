# R2 inner-law verification result

The inner-law operator checks pass, with one diagnostic explicitly unresolved.
The final Linux replay completed **4,366 independent checks**, covering 38
compiled kernel cases and eight retained flat 2D native operator fixtures. The
compiled kernel probe separately reports 128 native checks. The replay made
**zero native calls**, solved no implicit stage and advanced no trajectory.

The producing C++ images are from
`9f20b3053c303e8805f5183f668a8ec468d9656b`; the final Python analysis is
`ecfc04983d97e123760407c98d099f789ea208d4`. Both sources were clean and the
producing runtime fingerprints remained unchanged. Execution used CPU 50,
one rank and all six numerical thread controls set to one. The authoritative
report remains at
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r2/epoch-002/inner-reanalysis-004/results/validation.json`,
SHA256 `223be5d6ddfa79e07bd4460feab84c930f1650715490d6e34215cb75210a2986`.
The [compact machine result](inner-verification-result.json) binds the raw
fixtures, XML, binaries, imported sources, wrapper, log and runtime manifests.

| Check | Largest recorded component-scaled error (rounded) | Retained threshold |
|---|---:|---:|
| Compiled indexed stress | 5.55e-17 | 2e-12 |
| Compiled open-branch full stress action | 2.78e-17 | 3e-11 |
| Complete native alpha=0 weak loads, every row | 4.17e-13 | 2e-10 |
| Same-stress energy cancellation | 8.18e-14 | 3e-10 |
| Actual lifted inner traction samples | 1.11e-14 | 3e-10 |
| Mixed local-k/sampled-K/pressure/trace actions | 1.92e-9 | 8e-6 |
| Retained capped full weak + mass + constraints | 3.64e-14 | 3e-10 |
| Retained capped projected matching trace | 1.21e-15 | 3e-10 |

These are numerical implementation checks on the declared finite fixtures.
Their errors measure independent agreement, rather than small physical
residuals or zero interface jumps for these arbitrary verification fields.
The mixed-action threshold is the B component comparison; it does not replace
the separate C pressure-action physical-norm threshold of **1e-7**. That
separate pressure gate also passes and is hash-bound by the compact result.

The compiled checks include Linear and SaturatedR2 values, all nine gradient
directions, mixed local-k/auxiliary-viscosity dependencies, weak-divergence
production, full mean stress/traction, active and inactive branches, admissible
one-sided origins and invalid-input rejection. Inactive native fields and
unmasked operators agree exactly between the two selectors. The same selected
stress supplies momentum and actual production
`P=F:G-(2k/3)div(U)`. Scalar diffusivities, prescribed preprofiles, capped
normalization and the dependent EPS trace remain unchanged. Both pressures
retain their reduced convention; the outer molecular viscosity remains zero.

The legal cubic-depletion chi=0 fixture has rK=1, H=2 and lambda=1: it is
explicitly **uncapped**. The separately retained chi-active pressure fixture
exercises cap attenuation on 10 of 56 rays at those same physical heights,
with minimum lambda=0.5847870174412901. Its complete
alpha=300 operator is reconstructed from independent alpha=0 weak loads plus
the exact mass term and essential rows. Actual cap rays, pre/actual budgets,
and the baseline-plus-correction projection of lambda*Epre(H) pass. This is
not a cap-active chi=0 native test and is not a converged implicit stage.

The prescribed k~r² zero fixture retains positive sampled K and exact zero
local k on the ground. The actual covariance diagnostic evaluates 21 zero-k
points: SaturatedR2 reports minimum and maximum covariance **exactly zero**;
the matched Linear control reports the range
[-0.005239191889286443, -0.005204380946781486]. No outward central
k perturbation is taken there; compiled one-sided action checks cover the
admissible boundary. The fixture defines the zero trace before evaluation,
without a post-step repair or floor.

One statistic remains **unresolved at its original 3e-10 tolerance**: the
coefficient-derived maximum a/c in the saturated zero fixture differs by
2.278586980488665e-5 relative. At its native witness, r is about 4.85e-7,
k about 1.20e-13 and a/c about 9.17e10. The five other saturation statistics
retain their strict comparison. Independently checked exported point inputs
reproduce the actual selected-law quotient within 1.67e-15 relative. The
replay records primitive discrepancies and a directed-rounding, a-posteriori
same-point quotient bound; it does not certify all native basis errors or the
global maximum. The separate Linux 80-digit polynomial diagnosis is a finite
reference calculation, not an interval certificate. The machine result
explicitly sets `all_saturation_statistics_independently_certified=false`.
Nonzero-fixture statistic comparisons are unchanged.

Earlier failures remain retained: projected-versus-lifted state comparison,
invalid rK=H fixture, historical-oracle nested-exponential overflow and failed
NaN serialization, and the zero-fixture quotient comparison. The R2-specific
stable outer oracle changes no q, field or polynomial tensor. For q>8 its
omitted coefficient factor contributes less than 3e-64 in C and 4e-62 in C'
on the positive real axis; ordinary contraction roundoff is tested separately.
Historical M1/R1/M8 oracle sources remain unchanged.

This closes the focused inner implementation verification with the stated
statistic uncertainty. It establishes neither an admissible coupled stage
root nor physical outer-box coverage, polynomial-everywhere or fully discrete
positivity, curved MMS, long-time stability, scalability or terrain validation.
**M8 remains open.**
