# M8-R8 step selection and direction at the R7 endpoint

The reviewed source is `0ed7790450bb7e81f17ff1d0b577f12fc95f9cce`.
The saved iterate 25 is a nonlinear initial guess for alpha=400 at
time .0025 with the original h0 and startup history, eta=0 and zero
prescribed forcing. The 2,165-unknown outer6 p3/p1 spaces and quadrature,
SaturatedR2 with rho=.8, beta=.4, GammaScale=1, GammaUniform=0, original
outer response, viscosity and scalar diffusivities, projected epsilon
trace and capped source with theta=.2 and delta=.1 remain fixed.
The root, independent action and original-row feasibility gates stay
at 2e-12, 1e-7 and 5e-12.

## Bound experiment and allocation

`RANSWallR8Experiment=On` selects the fixed-state probe with
`RefinedR6`. The experiment reuses J25, the saved direction and its
action and the four native trial payloads on that ray. The loader
checks state, history, residual and matrix hashes, verifies the
retained trial states lie on the bound ray and compares the retained
action with the native J25 product. The base is already inside the
original box. True contact, the 5% scalar and continuous-wall reserves,
both complete original-row scans and the actual qualification rules
remain in use. Positive-value ground nondepletion is disabled for
the re-solved direction.

The hard local cap is 20 additional native trial states. Ordinary
damping can reach four fractions, from .0625 through .0078125.
Up to two smaller diagnostic fractions are available if all four
fail. The curvature selector shares the .0625 evaluation with
ordinary damping and can evaluate one independently selected fraction.
Its fallback reuses the completed ordinary results. Each attempted
radius is constructed once, with one informed adjustment reserved
for construction failure or a nonbinding radius. The informative
radius has at most five native trial slots. In-box partial outer
restoration cannot change the zero target. Reused states used for
source decomposition incur diagnostic cost without adding a new
trial state.

The predictor uses the retained .125 residual to fit the quadratic
correction to the entire residual vector. The held-out sample is
.0625, designated before scoring. Predicted values use the existing
block norms. The selector minimizes their maximum on a dyadic grid
from .0078125 to .0625. It halves the upper bound if held-out error
exceeds the larger of 2e-6 and half the actual held-out reduction
magnitude. The actual residual decrease gate is $10^{-4}f$ for
every candidate, with the actual coupled trace and original box
tested at each trial.

A same-ray reduction below .1% triggers a re-solve with radius four
times the accepted dimensionless displacement, bounded above by
.0625 and screened against twice the necessary algebraic radius.
The rowwise lower bound is a necessary condition. Full linearized
algebraic restoration and original-row checks still determine
whether a direction is usable.

## Local direction and finite-step behavior

The base reproduces the saved state, residual, constraints and history
byte for byte, including $E=2.926121942402$ and the positive ground
minimum $1.440972534\times10^{-5}$. The matrix hash is the saved
complete J25 hash. The independent action error is 2.049e-15 against
1e-7, using the saved R7 probe and the actual retained direction.
The scoped check supports the current action without a rebuild.

| Ray fraction | Linear-model $E$ | Native $E$ | Scaled vector remainder in k |
| ---: | ---: | ---: | ---: |
| 1 | 2.903696823 | 10.044679822 | 9.61548 |
| .5 | 2.909813711 | 3.513993471 | 1.80525 |
| .25 | 2.916702089 | 2.973295443 | .377523 |
| .125 | 2.921096958 | 2.931201796 | .0884604 |
| .0625 | 2.923530879 | 2.925883530 | .0215749 |
| .033203125 | 2.924725899 | 2.925384567 | .00602889 |

Ordinary damping accepts .0625 on its first added trial. Its actual
reduction per unit fraction is .0038146 against the model's .041457.
At .033203125 the actual and linear reductions per unit fraction are
.02221 and .04205. The k-vector remainder divided by fraction squared
is 5.66, 5.52 and 5.47 at .125, .0625 and .033203125. The values
support finite-step curvature with actual local descent. The recorded
profile and cap patterns are unchanged across the smaller ray trials.
The saturation-active volumes at base, .0625 and .125 are .056613568,
.056613569 and .056608747. These aggregate volumes do not identify
individual switched quadrature points.

At the held-out .0625 trial the predicted and actual maxima differ
by 8.04e-5. The k-vector model error is 7.98e-4, so agreement of the
two maxima is a weaker result than agreement of the complete residual
vectors. The selected .033203125 predicts $E=2.925416258$ and gives
$E=2.925384567$, with quadratic agreement 1.04491. The predictor
improves actual progress on the same ray and uses one more native
trial than the completed ordinary control.

| Candidate | Accepted fraction | Actual reduction in $E$ | Ground minimum | Trial cost |
| --- | ---: | ---: | ---: | --- |
| Ordinary continuation of retained ray | .0625 | .000238412 | 1.35542e-5 | 1 new evaluation, 1.849 s |
| Quadratic selector on retained ray | .033203125 | .000737375 | 1.39552e-5 | Held-out evaluation plus 1 new evaluation, 3.710 s, and .055 s model scoring |
| Re-solved smaller direction | .0625 | .000644146 | 1.35542e-5 | 5 evaluations plus 26.120 s direction construction |

The smaller radius .0166015625 is four times the best same-ray
displacement .004150390625. The necessary algebraic bound is
.0001874533. The constructed direction binds the new radius,
passes both original-row scans and has full linearized algebraic
residual below 1.4e-19 in the scaled blocks. It is a new constrained
direction with a different modeled residual, and its full actual
trial is 4.80909. Its first accepted trial is weaker than the
curvature-selected trial on the retained direction. The result
does not justify changing the radius floor in the pilot. No
second radius was needed. Seven additional native trial states
were evaluated in the whole fixed-state comparison.

## Component check and selected local qualification

The component audit uses the base and the existing .0625 and .125
states. It reconstructs the native residual within 2.50e-12 in
the scaled block norm. The following values are second divided
differences in the k block multiplied by .0625 squared.

| Term | Scaled magnitude |
| --- | ---: |
| Time mass | 1.09e-13 |
| Conservative advection | 3.48e-10 |
| k diffusion | 7.51e-11 |
| Inner production | 4.68e-7 |
| Inner destruction | .0226700630 |
| Outer production | 8.25e-9 |
| Outer destruction | 1.47e-15 |
| Complete constrained residual | .0226705117 |

The unchanged time-mass contribution is affine. The nonlinear
inner destruction load dominates the measured k curvature on
this ray. The result is a divided-difference observation and
does not identify a singularity or an exact termwise Jacobian
remainder. The component and branch diagnostics took 46.99 s
on reused states.

At the curvature-selected state, the momentum blocks are .11849118
and .22034576, the epsilon block is 5.09e-5, and the actual retained
trace differs from the projected target by .000845484 in scaled
coefficient units. The weak-continuity and gauge blocks are at
roundoff. The root threshold remains unmet. The governing gradient
clearance is 3.714e-8 at (-2.66666666667, 0, 2). The retained and
projected epsilon clearances at (1.75, 0, 2) are 8.107e-8 and
8.125e-8. Both pass the original box. The finite scalar inventory
and the continuous ground trace pass, with the minimum at
(1.56383558137, 0, 0). The measured covariance minima are
.00012075284 and .17687696020. Minimum lambda is one. The native
source defect E minus D minus C is about -1.3e-14 and same-state
energy weak-minus-direct work is -8.59e-13. This is a frozen trial
state with zero accepted physical steps.

## Selected short pilot

The same-ray selection gives 3.09 times the ordinary control's
residual reduction with roughly half its ground depletion. The
absolute improvement is still only .0252% of the base residual.
This supports a short pilot on residual progress per costly
current-direction solve. The pilot tests whether that gain
persists across fresh directions.

`RANSWallR8Experiment=Pilot` reproduces the first step with the
bound current J25, direction and two samples. After each accepted
state, the method constructs a fresh Jacobian and one radius-.125
direction with live fitting and validation samples. Up to eight
dyadic fractions from .125 and one model proposal are reachable,
with repeated samples cached. The first state already has its
two samples. After ten accepted updates or a completed-state
boundary at 300 solver seconds, the method publishes its final
report. The independent native wall cap is 450 seconds. The
original history remains fixed throughout.

The first selected state and residual reproduce byte for byte. Nine
new current matrices and nine new directions follow the bound first
step. The native process stops cooperatively at the completed
iterate 10 boundary and exits successfully after 347.10 s. The last
published current matrix belongs to iterate 9. The hard wall cap
and trial allocation were both available at the stop.

| Accepted update | $E$ | Ground minimum | Accepted fraction |
| ---: | ---: | ---: | ---: |
| 0 | 2.92612194 | 1.44097e-5 | — |
| 1 | 2.92538457 | 1.39552e-5 | .033203125 |
| 3 | 2.92362381 | 1.28203e-5 | .0478515625 |
| 5 | 2.91777683 | 9.95631e-6 | .125 |
| 7 | 2.91500137 | 8.51353e-6 | .125 |
| 9 | 2.91230420 | 7.15802e-6 | .125 |
| 10 | 2.91164151 | 6.73302e-6 | .0625 |

The pilot used 23 native trial evaluations and reached at most
three at any base. Eighteen trials passed the eligibility checks,
and ten were selected and saved as accepted updates. The selected
updates came from five model proposals, four live fitting samples
and one live held-out sample. The remaining five evaluations
failed actual residual decrease. All accepted states pass the
original box. Both original-row scans stay within 5e-12, with a
largest violation of 4.676e-12 across the fresh directions.

The residual falls .4949% across the pilot and the ground minimum
falls 53.27%. The ground ratios agree with $1-.95f$ within 1.06e-6.
Improved step selection helps the first fixed direction, while
the fresh-direction sequence remains slow and follows the
continuous-wall reserve. The cooperative limit is the observed
stop, with further admissible descent established and stationarity
unresolved. Assembly, factorization, direction construction and
model plus trial work took 87.16, 1.94, 155.05 and 46.97 s,
respectively.

The fixed comparison was produced at
`47d8e117e43d6f619b99c6d9e3f7401118cabf8c` and took 137.21
native seconds including its diagnostics. The pilot executable was
produced at `003714f33b21309295fdaf5484b6688a59aa1e08`. The first
launch preflight compared the local copy's pin with the producing
matrix path and stopped before native execution. The corrected
binding retains the identical J25 bytes and points to the original
payload. Both isolated Linux builds used 16 compiler workers. Native
solves used one rank on CPU 0 and one thread each for OpenMP,
OpenBLAS, MKL, BLIS, vecLib and NumExpr. Independent work used
CPU 49. The focused time, evolution and wall-stage native tests and
nine predictor binding tests passed on both builds.

## Independent pilot endpoint

| Separate outcome | Terminal observation |
| --- | --- |
| Coupled root | Independent maximum 2.91164151002389 against 2e-12 |
| Actual retained and projected trace | Scaled matching discrepancy .000424296895 across 16 rows |
| Original box | Passed. Gradient clearance 1.721e-8 and retained and projected epsilon clearances 3.756e-8 and 3.765e-8 |
| Scalar positivity | Passed on the original finite inventory, with zero independently unresolved signs |
| Continuous ground positivity | Minimum 6.733016805e-6 at (1.56384103496, 0, 0), confirmed by 80- and 110-digit reconstruction |
| Covariance and lambda | Minima .00011926452 and .17687700000. Minimum lambda one |
| Direction accuracy and model agreement | Original-row gate passed. The verified J25 action supports the initial curvature comparison, with fresh live models on later states |
| Physical history | Startup slots verified against the native final report. Zero accepted physical steps |
| Stop | Cooperative completed-state limit with final report published |

The momentum blocks are .11848725 and .22024526, the k block is
2.91164151 and the epsilon block is 5.125e-5. The weak-continuity
and gauge blocks are at or below 2.37e-15. The independent maximum
scaled block difference is 2.99e-12. The maximum matching trace
coefficient difference is 3.082e-5 before epsilon scaling. The
governing gradient witness is (-2.66666666667, 0, 2) in element 31.
Both terminal epsilon witnesses are at (1.75, 0, 2) in element 34.

The independent source-budget defect E minus D minus C is -1.65e-15
and same-state energy weak-minus-direct work is 1.06e-13. The native
and independent energy diagnostics agree to 1.07e-12. The normal
and tangential traction jump RMS values are 3.83e-4 and 2.06e-4.
The k-flux jump RMS is 6.53e-4. The accounting identities concern
the saved non-root.

The [decision ledger](cls-epawf-revised-evidence/m8-r8/decision-ledger.json)
pins the replay matrix, retained trials, current action check, local
and pilot runtimes, endpoint and independent audit. The full
payloads remain under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r8`.

## Next decision

At the saved in-box R8 endpoint, assemble its current matrix and
compare one coupled direction model that includes nonlinear inner
destruction with the existing linear construction. Match the wall
reserve, box tests and work allocation. Judge the actual residual
gain against ground depletion. Step selection already improves the
retained-ray result, and the pilot shows that this gain alone leaves
slow convergence and substantial depletion. The measured inner
destruction curvature gives a concrete term to target in the
direction model while retaining the same physical equations.

M8 remains open. Curved-wall MMS, active-cap curved or variable-height
coupling, full-element and fully discrete positivity, refinement,
long-time robustness, scalability and physical floor/hill validation
remain outstanding.
