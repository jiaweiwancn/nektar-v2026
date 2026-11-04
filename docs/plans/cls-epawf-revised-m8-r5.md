# M8-R5 joint residual and outer-region restoration

The exact target retains the SaturatedR2 law with saturation parameter .8,
the original outer response with beta=.4, GammaScale=1 and GammaUniform=0,
the original molecular-viscosity conventions, scalar diffusion, projected
epsilon trace and capped source with theta=.2 and delta=.1. It uses the
2,165-unknown outer6 p3/p1 representation and its existing quadrature. The
3,000-unknown guard remains in force. The stage is alpha=400 at time .0025
with the original initialized h0, eta=0 and zero prescribed forcing. The
root and action thresholds remain 2e-12 and 1e-7, respectively.

The reviewed branch was `d3c236ad027caf2a69b1d84db942f1943109b9c4`.
Linux built `abc370f66414ecfc0a00ce8ab8628cd88a75ec1b` in an isolated
checkout with 16 compiler workers and the existing MPI Release libraries.
The solver used one rank on CPU 0, with OpenMP, OpenBLAS, MKL, BLIS, vecLib
and NumExpr thread counts set to one. The canonical Linux source checkout
and the historical producing builds were preserved.

## Terminal replay and box witnesses

The current executable reproduced the R4 terminal iterate 12 state,
original h0, exact coefficient residual and constraint payload byte for
byte. The complete residual maximum is 4.995492353429328 in k. The finite
scalar inventory passes. Reconstruction of all five ground cubics at
80 and 110 decimal digits resolves the minimum as positive,
.0016375755634 at x=1.56414545965. The measured covariance minima are
.0004378988654 inner and .1770776315 outer. The observed minimum lambda is
one. The retained R4 matrix belongs to iterate 11. The matched local
experiment assembled and factored a new current matrix at iterate 12.

| Quantity | Location and observation | Actual value | Fixed bound | Positive normalized excess |
| --- | --- | ---: | ---: | ---: |
| (G_{xz}/g_0-G_{0,xz}/g_0) | (1.66666666667, 0, 2), element 34, outer-side matching transfer | -.0035647171382 | -.001 | 2.5647171382 |
| Retained epsilon/reference | (1.75, 0, 2), element 34, outer-side dense reconstruction | .9744018120862 | .98 | .2799093957 |
| Projected epsilon target/reference | Same matching-interface coordinate and outer element | .9744022579708 | .98 | .2798871015 |

Both epsilon observations violate the lower bound. The retained polynomial
gives the lower of the two. Outer k and sampled K remain within their
unchanged intervals. The projected and retained ratios are distinct because
the actual trace constraint at this non-root is nonzero.

## Direction model and nonlinear gate

`RANSWallGlobalization=JointRestorationR5` opts in to the numerical policy.
`RANSWallStageTest=DirectionExperimentR5` compares the two directions at a
single bound state. `ExactTargetR3` with that globalization and an explicit
`RANSWallR5Updates` value selects the continuation. The older selectors
keep their original paths.

For each measured quantity (u), use lower and upper excesses
((l-u)/w) and ((u-r)/w), where (l,r) are its original limits and
(w) is the permitted half-width. The half-widths are .001 for gradient
departures, .01 for outer k and sampled K ratios, and .02 for epsilon
ratios. (V) is the largest positive excess, with zero for a passing
family. The report also records each of the four family maxima and the sum
of positive excesses across the inventory. That sum retains duplicates
across observation types. There are 54,385 quantity records in the base
inventory. The flat 2D velocity has four represented gradient components.
The other tensor components are identically zero for the frozen reference
and spaces.

The rows use the original audit locations, one-sided element gradients,
reference scales, and both actual and projected epsilon values where the
audit checks them. Outer-bound locations are fixed on this fixture.
Normalization rays are open and inner. Their moving abscissae and weights
still affect the projected-trace derivative through the verified current
constraint action. For a trace coefficient whose residual is the retained
value minus its projected target, the projected action is its retained
step minus that row of (Jd). Nonlinear trial evaluation rebuilds the
scalar inventory, external lift, trace projection and normalization.

The direction solves the R4 residual-space least squares
[
  d=d_N+J^{-1}Bz,qquad d_N=-J^{-1}R,qquad
  min_z 	frac12|z|_2^2 .
]
Each free differential block of (B) uses the original consistent-mass
Cholesky factor multiplied by its physical rate scale and square root of
volume. The algebraic rows of (B) are zero. This gives
(C+J_Cd=0), including the nonzero base trace defect, weak continuity,
gauge and external constraints. A twice reorthogonalized active-set QR
solves the cut problem with dual releases. The direction retains the R4
dimensionless infinity radius 2, the 5% linear scalar and continuous
ground reserve, and the nondepletion rule at near-active ground extrema.
The actual nonlinear checks use fresh ground extrema and the original
finite scalar inventory.

The first model requests zero excess in every family. If the bounded
full exchange fails outside the box, the partial model can reduce each
positive family maximum by 25% and retain the original bounds in passing
families. These ceilings apply only to the direction model. Neither adds
a PDE term, wall condition or acceptance slack. The local gate requires
at least 1% actual decrease in both the complete residual and (V), with
no worsening of any family maximum. The continuation requires residual
decrease and a reduction of (V) by at least (10^{-4}fV) outside the box
at trial fraction (f), with every family maximum nonincreasing. Once
inside, it requires every family to remain inside. No temporary residual
growth is accepted.

## Matched fixed-state comparison

| Direction and trial | Predicted (E) | Actual (E) | Actual gradient excess | Actual epsilon excess | Ground minimum k |
| --- | ---: | ---: | ---: | ---: | ---: |
| Replayed base | — | 4.99549235 | 2.56471714 | .27990940 | .00163758 |
| R4 control, full step | 3.10670339 | 3.92188670 | 2.65236996 | .32392410 | .00151711 |
| R5 partial restoration, full step | 3.10910963 | 3.92879913 | 1.92353785 | .20993205 | .00151804 |

Both directions pass the finite scalar checks. The partial trial lowers
(E) by 21.35% and both positive violation maxima by 25%. Its modeled and
actual box maxima agree within 1e-12, while the model overpredicts the
residual reduction. Its measured covariance minima are .0003039224 and
.1770604981, lambda is one, and the recorded profile branches are
unchanged. It remains outside the box. The sum of positive excesses
changes from 1,070.175 at the base to 1,113.385 for R4 and 1,054.636 for
R5. The exact complete block norms are in the decision ledger.

Full linear restoration exhausted 96 cut scans with a remaining normalized
cut violation of 4.65e-11 against the direction-model threshold 5e-12.
The last reported cut was an outer native-reconstruction gradient row
at (-2, 0, 2). This is a numerical exchange limit. The partial solve
converged in 24 scans with 17 releases. Its active rows include the
near-active ground extremum, a continuous reserve row, three outer
gradient rows and the projected epsilon target at (1.75, 0, 2). The
scaled projected-row multiplier is .14075638. The infinity step is
1.33942, below radius 2. The ground and reserve rows are numerical
safeguards, and neither trial contacts a scalar evaluation boundary. Stationarity, active equality error and
complementarity were 1.72e-15, 1.70e-16 and 1.14e-15. The maximum
linearized algebraic residual was 2.26e-14. The actual trial still has a
.000421899 scaled constraint block, whose maximum occurs on a matching epsilon
row. The useful joint descent was obtained
with the ground safeguards enabled, so there was no matched removal of
the inward rule at the base.

The independent R3 value oracle reused the single bound current matrix for
both directions. At the two complex-step sizes its maximum physical
action errors were 1.21e-14 for the R4 control and 1.22e-14 for R5 partial
restoration. Independent normalized derivatives of the controlling box
witnesses, including the projected epsilon value, agreed within
5.43e-16. The inherited guards for exact cap knots and moving points
outside their selected element did not fire at assembly. These checks
qualify the tested directions at the replayed state.

The measured full action assembly took 9.158 s, factorization .188 s and
the three direction calculations 19.522 s combined. The pilot budget was
set before launch at 19 updates, 32 trials per update and 900 native
seconds. The estimate was 41.244 s per update including six trials at the
local mean cost, plus 90 s for setup and terminal checks. The independent
wall cap covers variation in the actual trial and exchange costs.

## Single exact-target continuation

The first and only accepted update reproduced the partial-restoration
trial byte for byte. The initial current matrix also reproduced the
local comparison matrix. The solver then assembled and factored a second
current matrix bound to the accepted state. Both full and partial
direction solves exhausted their bounded exchanges at that state. The
last gradient-cut violations were 3.57e-11 for full restoration at
(-3, 0, 2.33333333333) and 9.02e-12 for partial restoration at
(-1.58333333333, 0, 2). The reported active restriction was an outer
gradient cut in both cases. No second nonlinear trial was evaluated.
The run used 98.16 native seconds. The stop was the numerical exchange
limit, with 18 permitted updates and most of the wall budget unused.

| Separate outcome | Saved terminal observation |
| --- | --- |
| Exact target convergence | Failed. Native and independently reconstructed maxima are both 3.92879913247 in k against 2e-12 |
| Finite scalar admissibility | Passed on the original declared inventory |
| Continuous ground safeguard | Passed. Minimum k .00151803603 |
| Fixed constitutive region | Failed. Gradient departure .00292353785 and minimum epsilon ratio .97580135906 |
| Measured covariance | Positive minima .00030392243 inner and .17706049812 outer |
| Lambda qualification | Passed. Minimum one |
| Physical history advancement | Zero steps. The stage remained an unqualified numerical iterate |

The native weak-continuity blocks are 8.78e-16 and 8.34e-16 and the gauge
block is 1.08e-17. The separate constraint block is 4.22e-4. The near-zero
linearized algebraic defect did not carry over to the finite nonlinear
trial. The existing certificate still excludes the accepted state.
Measured covariance does not establish incremental ellipticity there.
The small violations reported by the direction solver concern the linear
search cuts at its last exchange, whereas the actual saved state has
(V=1.92353785). The continuation therefore leaves local compatibility with the fixed
physical region unresolved.

The maximum native-to-independent scaled residual difference is 3.15e-12.
The independent physical inventory agrees with the native applicability
and covariance checks. The independent source defect E minus D minus C
is 7.68e-17. Same-state weak energy work minus direct work is 5.42e-14
and the native-to-independent energy diagnostic agrees to 4.96e-13.
The normal and tangential total-traction jump RMS values are 3.91e-4 and
3.66e-4, and the k-flux jump RMS is .0011313. The accounting comparisons
are at the same retained non-root.

The three focused native time, evolution and wall-stage regressions and
nine predictor binding tests passed. The durable collector was also
checked against the R5 progress stream and the historical R4 stream.

The [decision ledger](cls-epawf-revised-evidence/m8-r5/decision-ledger.json)
binds the source, executable, XML, initialized history, compared
directions, both pilot matrices and the terminal audit. The Linux payload
root is
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r5/epoch-001/cases-001`.
It retains the full comparison, ground extrema, trial state and residual,
independent audits, and durable stage progress.

## Next decision

At the saved R5 terminal state, condition the nearly repeated gradient-cut
rows and independently measure their attainable feasibility and KKT
accuracy before launching another stage continuation. Exact duplicate
coalescing and a higher-precision verification of the projected active
solve can test the cause of the exchange limit while retaining every
nonlinear scalar, box and root gate.

M8 remains open. Curved-wall MMS, active-cap curved or variable-height
coupling, full-element and fully discrete positivity, refinement,
long-time robustness, scalability and physical floor/hill validation
remain subsequent gates.
