# Final bounded controls: both stop before the first converged stage

The [controlled ledger](controls-ledger.json) and complete compact
[dt=.005](controls-dt005-iterations.json) / [dt=.0025](controls-dt0025-iterations.json)
iteration records preserve both unsuccessful final allocations. Independent
analysis passes 12,632 norm, constraint, linear-product, scalar and
offending-point checks. It launches no solver or retry. Complete native
vectors, histories and Jacobians remain in the source-bound Linux traces;
their exact paths and SHA256 hashes are retained.

| Control | Native time | Newton iterations | Internal accepted descents | Residual exceptions | Valid non-descents | Completed steps |
|---|---:|---:|---:|---:|---:|---:|
| dt=.005 | 274.99 s | 9 | 8 | 95 | 0 | 0 |
| dt=.0025 | 275.79 s | 12 | 11 | 117 | 0 | 0 |

Both have one accepted physical state: their common initial state at t=0.
The initial coefficients and native startup-history values exactly match the
refined zero-step preflight. The dt-specific IMEX2 history labels are correct,
and the first implicit startup stages have alpha `200`/`400`, at times
`.005`/`.0025`. Neither startup stage converges. No physical guard stop,
completed time step, timestep-refinement result or physical passage is claimed.

The only saved physical-range event is the non-stopping
`outer_gradient_box_exit` at t=0. All attempted continuation is explicitly
unqualified. The [focused guard identity checks](controls-focused-guard-identities.json)
also confirm exact event-state/history and complete native-history identities
against the corresponding unguarded preflights. Strict `R1LocalV1` stops the
refined initial state solely for its local-box exit. The original coarse
`R1DiagnosticV1` zero-step check stops for its already negative covariance.
Neither focused check appends an accepted state or audit epoch.

## What rejects the final nonlinear searches

Every one of the 212 residual exceptions is independently located on a
matching-trace **GLL tangential ray**, during complete-ray capped-source
normalization. These are near-wall normal-quadrature points along a ray
anchored at the matching trace. They are neither negative k at H nor the
volume-source Gauss ray locations. Their transported k is negative while
the physical sampled K at rK remains about `.51053`.

The final rejected searches each evaluate all 21 unchanged fractions through
`2^-20`. At that smallest fraction, independent k is approximately
`−1.64477336e-7` and `−7.45760685e-8`, against arithmetic scales about
`9.1e-15`. The current base iterates at those same coordinates have positive
k about `1.93e-7` and `4.89e-9`; even the smallest permitted trial leaves
the finite normalization domain. No trial yields an evaluable residual that
merely fails the strict descent comparison in these controls.

Internal Newton iterates already have negative source-quadrature inner
covariance at iteration 1. Their native GLL k first becomes negative at
iterations 5 and 7, independently confirmed, while the volume Gauss source
k remains positive. Every record includes the actual full G/S, k, sampled K,
k/K, nuA and rR at the exported offending inner points, with formal algebra
explicitly labelled when k is negative. The same-point initial tensor is
retained to distinguish local depletion from the original coarse
representation's excessive strain.

These internally accepted *Newton descents* are not accepted physical
stages. The finite physical guard is reached only after a converged nonlinear
solve, which never occurs here. Native source-ray summaries retain lambda
one and zero cap-attenuation volume; this is not evidence that cap attenuation
caused the stop. A source-rate limiter does not enforce tensor realizability
or polynomial-everywhere scalar positivity.

The single initial Jacobian is reused through ages 8 and 11. The final scaled
k residuals remain about `26.83` and `11.33`; divergence/pressure and gauge
rows are near roundoff. All saved linear products and residual norms are
checked against the exported cached matrix. No fresh-J comparison or retry
was made for these final allocations, so refresh-policy causality is not
proved by them. The separate original-case controlled comparison remains
distinct evidence.

## Provenance and limits

The producing source is `360cbb7ac728b575cdf9af21e6d0963a1d88acb9`; the
executable SHA256 is
`aa48843d03d22b7089b2b08f56d777640a1419a3e59296850523d3130206e336`.
The [runner](controls-runner.json), [plan](controls-plan.json),
[analysis runtime summary](controls-analysis-runtime.json) and
[curation index](controls-curation-index.json) bind actual commands, XML,
source, numerical modules, histories and results. Both producing and analysis
before/after fingerprints match, retaining 3,921 inputs and 100 runtime
images. Analysis uses CPU 49 and one numerical thread in all six recorded
thread controls. These are matrix slots 5 and 6; no additional evolution is
permitted or performed.

The frozen backend allows inner k=0 but evaluates covariance only for k>0.
For positive sampled K, the original linear inner tensor is still defined
at k=0, so this is not a complete exact-zero covariance guard. No exact-zero
k sample occurs in the independently scanned stored-iterate native GLL sets;
negative internal samples do occur. The source is not changed after these
runs, and no clipping, field floor or relaxed tolerance is used.

The common initial state is finite-admissible but outside the frozen outer
region; the failed implicit attempts are unqualified diagnostics. No
completed unforced refinement/long-time result, continuum instability proof,
new coupled inner-law verification or physical terrain validation follows
from these controls. M8 remains open.
