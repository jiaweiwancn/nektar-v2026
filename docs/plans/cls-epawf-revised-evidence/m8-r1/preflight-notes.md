# M8-R1 representation preflight

Exactly four zero-step initializations compare the actual failed flow mesh
with one frozen radial refinement, at two source quadrature levels. All four
independent implementation audits pass, totaling 532 recorded checks.
The refined initial state has positive covariance and positive transported
scalars on every declared finite set. It remains outside the already frozen
outer gradient-entry radius `.001` region. No trajectory is run or qualified
by this evidence, and no new inner constitutive completion is adopted.

The [summary](preflight-summary.json) carries exact values, witnesses,
state/history/result hashes and quadrature comparisons. The
[completion runner](preflight-runner.json), [runtime summary](preflight-runtime-summary.json)
and [postprocessing receipt](preflight-postprocess-receipt.json) bind the
actual solver, reviewed `ab5c3c6fee19761e8d88692e8e6c331e5889234e` source,
imported modules, child environment and commands. The pinned historical
corrected solver has SHA256
`71897ba6d378acf35f84ae832920ed17cc3c66db8b9f193211e46c5ab043f1e4`.
All 1,674 checked C/C++ inputs agree between the reviewed archive and the
actual build source. The full before/after manifests retain 3,920 inputs and
100 runtime images, with no source, module, library, tool or environment change.

## Frozen representations and cost

Both use the physical domain `[-4,2] × [0,4]`, tangential knots
`[-4,-3,-2,-1,1,2]`, p3 velocity/scalars, p1 pressure and GLL8 points.
The original radial mesh is `[0,2,4]`. The sole refined mesh is
`[0,.01,.1,.25,.5,1,2,8/3,10/3,4]`. H, roughness, physical sampling
distances, transition, source parameters, beta `.4`, GammaScale `1`, theta `.2`,
delta `.1`, eta `0`, dt, tolerances and solver policy are retained.

| Representation | Cells | Full unknowns | Single dense array | Initialization, base/high |
|---|---:|---:|---:|---:|
| Original | 10 | 425 | 1.445 MB | 8.0 / 14.0 s |
| Refined | 45 | 1,571 | 19.744 MB | 10.0 / 20.0 s |

The unchanged dense guard is 3,000 unknowns. Relative to 425 unknowns,
the refined mesh has an indicative `16.63×` volume finite-difference
Jacobian assembly cost and `50.51×` dense factor work. These estimates apply
to a future nonlinear stage; the zero-step timings do not measure that cost.
Execution is sequential on CPU 49 with one numerical worker and all six
recorded numerical thread variables set to one.

## Initial physical and representation observations

Base source quadrature is tangential/normal `8/16`; high quadrature is
`12/24`. Space and GLL points remain fixed between quadrature levels.

| Case | Minimum inner covariance | Minimum outer covariance | Maximum outer normalized gradient departure | Outer epsilon/epsilon0 |
|---|---:|---:|---:|---|
| [Original/base](preflight-original-base.json) | −.09773632 | .17433082 | .01833233 | [.98325247, 1.00408123] |
| [Original/high](preflight-original-high.json) | −.09774725 | .17432242 | .01833233 | [.98325247, 1.00408123] |
| [Refined/base](preflight-resolved-base.json) | .07729332 | .17679689 | .00341234 | [.98325247, 1.00006559] |
| [Refined/high](preflight-resolved-high.json) | .07729332 | .17678911 | .00341234 | [.98325247, 1.00006588] |

The original/base initial coefficient vector is bitwise identical to the
archived failed case's first accepted t=0 vector, and the physics signature
is exactly preserved. Raising quadrature alone retains its negative inner
covariance. The maximum finite inner ratio
`r_R = 3 nu_A lambda_max(G+G^T)/(2k)` decreases from about `1.2872` to
`.7729` after the single spatial refinement; local k remains undepleted.
The maximum relative inner AWF gradient error decreases from `1.7385` to
`.6539`. Positive refined covariance is a finite initial-state observation,
not a polynomial-everywhere or evolving-state theorem.

Every case has zero limiter-active and cap-active volume and minimum lambda
one in the actual source diagnostics. The independent source normalization,
redistribution, production/stress exchange, complete traction, k flux and
dependent epsilon projection agree with the native diagnostics. The refined
trace projection error divided by E0 is below `5.32e-15`. Actual chi-zero trace
epsilon ratios span approximately `[.98345285,.98465195]`; a ±1% epsilon
certificate cannot be presumed to cover this baseline trace.

The worst refined outer gradient departure is one-sided at `r=8/3`, with
normalized `Gxz=1.0034123386947804`. Its native source Gauss sets also reach
`.00291812`, above `.001`. The frozen screening region is therefore not
covered even though its expanded epsilon bounds `[.98,1.02]` include these
samples. No mesh, state box, physical input or acceptance tolerance was
changed after observing this result. A later parent-authorized control must
be explicitly diagnostic and unqualified from t=0.

Finite sets include source and mass quadrature, matching and sensor/transfer
locations, native GLL reconstruction, outer dense `2p+3` nodes per axis,
and additional inner dense `2p+7` nodes per axis. One-sided element
gradients are retained. These checks do not establish a continuum invariant
region, fully discrete turbulence admissibility or physical terrain accuracy.

## Raw residual and quadrature sensitivity

The independent complete unanchored weak loads agree with the native loads
to `3.56e-14` in the refined/base case. Pressure/divergence and common-gauge
residuals are near roundoff, while differential residuals remain substantial.
The refined/base physical block norms are approximately
`[9.86323, .025532, 62.59848, .111044]` for velocity/velocity/k/epsilon.
This is an unforced initial state intended to evolve, not a stationary
reference solve. Mechanical-pressure continuity is not imposed.

Base and high initial coefficients are bitwise identical within each
representation. Increasing only the independent normal integration order
from 16 to 24 changes complete weak loads by less than `9.44e-15`.
Increasing tangential/outer quadrature as well changes native weak loads by
`7.06e-5` on the original mesh and `6.14e-5` on the refined mesh; the refined
k-block norm changes by `.00215343` on a norm near `62.6`.
The refined physical mass-dual norm of that k-load difference is `.00834036`.
Two quadrature levels provide a sensitivity observation, not an order or
convergence claim.

The [initial attempt](preflight-initial-attempt.json) preserves a tuple/list
serialization mismatch in the plan's time-integration metadata. Its
successful original/base initialization was retained; only its independent
audit was repeated after the metadata-only fix, followed by the three
missing initializations. Exactly four native initializations ran in total.
The [curation index](preflight-curation-index.json) pins all 18 byte-for-byte
copies; large native fields, full histories and runtime manifests remain
in the Linux run directories. M8 remains open, with genuine curved MMS,
active-cap curved/variable-height coupling, full turbulence admissibility,
completed unforced refinement/long-time studies, scalable solves and paired
developed terrain/reserved-score validation still separate gates.

## Prepared diagnostic controls

The source case generator's `--family trajectory` mode prepares only the
approved resolved/base controls: dt `.005` with 12 steps and dt `.0025` with
24 steps, both through at most T `.06`. It writes no body/scalar forcing,
retains the original physical/function strings and solver tolerances, and
sets restart/checkpoint to `None` so each dt initializes its own valid IMEX
history. Static XML previews pass these checks; evolution is parent-owned
and is not part of the four preflight results.

The exact opt-in solverinfo values are:

```text
RANSWallPhysicalRangeGuard = R1DiagnosticV1
RANSWallSolveDiagnostics = Failure
RANSWallSolveDiagnosticReplay = Off
RANSWallSolveDiagnosticFile = case.nonlinear.json
```

The guard records the first frozen-local-region exit at t=0 and labels any
permitted continuation diagnostic. Its separate finite physical gates stop
the first declared scalar/K, covariance or lambda violation before history
advance. The prepared plans introduce no numeric guard parameters or new
model law. Root freezes the generator, pins the new guard build and allocates
the two sequential evolution calls; the preparation function invokes no solver.
