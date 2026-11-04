# M8-R4 coupled feasible directions at the scalar boundary

The local experiment found a genuinely coupled feasible descent direction for
the unchanged SaturatedR2 alpha400 stage. Its accepted quarter step reduced
the scaled residual from 17.3238 to 14.0961. The safeguarded ordinary Newton
comparison reached 16.1365, and R3 had stalled at 16.0695. M8 remains open.

The reviewed baseline was `797ca6629299e8433c7dd524e958d8b7f11002a3`.
The exact problem uses the original initialized h0, time .0025, alpha400,
eta0 and zero prescribed forcing. The outer6 mesh has 2,165 unknowns,
p3/p1 spaces and the retained quadrature. The 3,000-unknown guard,
SaturatedR2 stress with rho=.8, regularized outer law with beta=.4,
GammaScale=1 and GammaUniform=0, viscosities, scalar diffusivities,
transition, profiles and projected dependent trace retain their R3
definitions. The capped destruction parameters are theta=.2 and delta=.1.
The action threshold is 1e-7 and the coupled residual threshold is 2e-12.

## Recovery and fixed outer region

I evaluated the exact binary64 coefficients in 80- and 110-digit Decimal
arithmetic. For each of the five ground elements I used the signed native
tensor map to obtain the four ground modes, converted them to a cubic in
the element coordinate on [-1,1], and evaluated both endpoints, every real
interior stationary point and the original limiting GLL location on its
element. The precision comparison agreed below 1e-65 in these values.
This check covers the ground trace in x.

| Saved R3 iterate | Native minimum inner k | Wall minimum | x at wall minimum | Maximum outer gradient departure |
| ---: | ---: | ---: | ---: | ---: |
| 0 | .344686 | .344630 | 1.613666 | .000532248 |
| 1 | .00172343 | .00152346 | 1.614391 | .00189601 |
| 2 | 8.61715e-6 | -6.82274e-5 | 1.610737 | .00196552 |
| 3 | 4.30857e-8 | -7.63215e-5 | 1.610718 | .00196588 |
| 4 | 2.15429e-10 | -7.63620e-5 | 1.610718 | .00196589 |
| 5 | 1.07714e-12 | -7.63622e-5 | 1.610718 | .00196589 |
| 6 | 5.35357e-15 | -7.63622e-5 | 1.610718 | .00196589 |

For iterate 4 on [1,2], the signed coefficients in the Modified A basis
are `(.03444218891017109, .5105080130598861, -1.005332593665448, -.8191849353388675)`.
The endpoints give the first two coefficients. The stationary values are
`.1265613695716572` at `1.1847430898946725` and
`-7.63619554992865e-5` at `1.6107177089440379`. The original GLL value
reconstructed at `1.6046496089512394` is positive at `2.154286817620477e-10`.
On the adjacent element [-1,1] the coefficients are
`(.07408665923481723, .03444218891017109, .0802254175906373, .06919127829708979)`.
Its interior stationary values are `.0656963544470869` and
`.07623677367849357`. The full endpoints and station locations for every
trace are in the retained wall check.

Checking backwards from iterate 4 selects iterate 1. Its limiting-element
coefficients are
`(.037481109324371154, .5105080130598861, -.9978092176697364, -.8319412142631637)`.
Its minimum is positive at every checked ground element. R3's finite scalar
inventory classifications retain their original scope, including the
sign uncertainty reported for the terminal native sample. The R4 seed
retains all original coefficients and its nonzero trace residual. It is
an internal numerical guess with the original h0 in the history slot.

Iterate 1 is the first observed exit from the fixed outer box. At
`(1.6046496089512394,0,2)`, element 34 on the outer side, the normalized
`G[x,z]` component departs from reference by `-.0018960121`, beyond the
.001 bound. The other scalar ratios remain in their fixed intervals in
the saved records. Both covariance minima at iterate 1 are positive
(.0003363705 and .1770926713). The selected seed is an unqualified algebraic
research state.

## Coupled direction

`RANSWallGlobalization=CoupledR4` opts in alongside `ScalarFeasibleR3` and
`VerifiedR3Action`. The ordinary Newton step is $d_N=-J^{-1}R$. For each
free differential block $b$, put its scaled Cholesky factor
$(rate_b\sqrt{V_b})L_b$ into $B$, where $M_b=L_bL_b^T$ is the existing
free consistent mass. Its algebraic rows are zero. Search in the form
$d=d_N+J^{-1}Bz$. Thus $R+Jd=Bz$, and all weak continuity, gauge,
external and actual projected dependent-trace equations satisfy
$C+J_Cd=0$ to linear solve accuracy. The objective is
$\lVert z\rVert_2^2/2$, which equals the sum of squared physical
scaled differential block residuals divided by two. There is no
regularization term. The dimensionless infinity radius is 2, using
$u_*,u_*,K_0,E_0,u_*^2,u_*^2,u_*/H$ by coefficient block.

The frozen current scalar rows and continuous ground trace must retain at
least 5% of their base values in the linear model. The original limiting GLL
row and all current ground extrema at or below $.005K_0$ must have
nonnegative directional changes, which permits inward motion. The value row
at a stationary point supplies the derivative of that extremal value.
Exchange uses the most violated current cut, including extrema of the
modeled reserve polynomial. The projected dependent-trace inequalities use
the actual constraint action to recover the projected target derivative.
This numerical safeguard operates in the search. Each nonlinear trial
rebuilds moving rays, external lifting and the dependent trace through the
existing scalar-domain routine and checks fresh ground extrema.

The active solve works in residual coordinates using the current pivoted
factorization of J and twice reorthogonalized QR of the small active normal
set. It tracks primal feasibility, nonnegative dual multipliers,
stationarity and complementarity. A dual step keeps the entering cut
through blocking releases. This also permits return to ordinary Newton
when it meets the trust and scalar conditions. The first local version
failed on a dependent active row. A rank-aware exchange alone still
exhausted the 40-pass bound. Both local failures are preserved. The
completed version uses the persistent entering cut, which resolved the
dependency cycling. The same C++ direction calculation is used in the
local check and the pilot.

## Local comparison

The base and rebuilt complete matrix match the retained iterate-1 state,
residual and Jacobian bytes. The independent R3 complex action machinery
checked the actual R4 direction at steps 1e-20 and 1e-30. The maximum error
was 1.73e-14 against 1e-7. The R4 assembly checks reject a rounded cap-knot
equality or a moving row outside the selected element. Neither guard fired
during local qualification. The trial reports record branch changes and
the actual merit rather than extrapolating through them.

| Step | Fraction | Native maximum scaled residual | Ground/finite scalar checks |
| --- | ---: | ---: | --- |
| Base | — | 17.32380644 | Pass |
| Safeguarded Newton | .06598103 | 16.13652749 | Pass |
| R4 full | 1 | 34.50176079 | Pass, residual increase |
| R4 half | .5 | 27.37237927 | Pass, residual increase |
| R4 quarter | .25 | 14.09613549 | Pass, accepted |

The ordinary direction changes k at the limiting GLL row by -.02476968.
R4 changes it by +.001001422. The dimensionless R4 full step is 1.84393.
Its linear model maximum is 3.66713 and its algebraic error is below
1.8e-16. It used 16 cut scans and 12 releases, four of them dependent
exchanges. The active rows were the native point x=1 and the continuous
reserve rows x=-2.1503754 and x=1.5632164. Their multipliers in the
declared constraint scaling are 18.2513, 42.8514 and 222.4912.
Stationarity, active equality error and complementarity were
5.74e-16, 4.44e-16 and 1.55e-15. At fractions 1, .5 and .25 the linear model predicts 3.66713, 9.43881
and 13.32967. All three nonlinear trials change the recorded profile
branch pattern. Their minimum b stays above .824 and lambda remains one.
The feasible line search evaluates the original complete maximum norm. At the accepted trial the projected-trace
constraint block is .07452, even though the undamped linear model restores
it. The outer departure also increases to .00278858. Covariance remains
positive on the measured inventory (.000386888 inner, .176995556 outer).

The measured direction cost was 2.109 s and full action assembly 9.164 s.
The 12-update budget projection, with up to 32 trials per update, was
422.2 s plus a 90 s setup and audit allowance. The independently verified
quarter-step reduction passed the predeclared 1% improvement criterion
against the base, safeguarded Newton and R3 terminal maxima. The pilot
therefore uses a 600-second native wall cap from the same saved seed.

## Exact-target pilot

The producing Linux revision is `a61a607d424d003323bdc69a80e5543daf0896fb`.
The build used the installed MPI Release libraries and 16 compiler workers.
The pilot used one rank on CPU 0, with OMP, OpenBLAS, MKL, BLIS, vecLib and
NumExpr thread counts all set to one. The producing checkout, runtime, XML,
h0 and all durable payloads are pinned in the decision ledger. The three
focused native time, evolution and wall-stage regressions passed. The
initial exact residual and rebuilt action matrix reproduce R3 at the
selected seed. The regression and replay evidence concerns the retained
equilibrium operator and solver entry points. The pilot concerns the
nonlinear algebraic stage.

| Completed update | Native maximum scaled residual | Ground minimum k | Outer gradient departure | Accepted fraction |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 17.32380644 | .00152346 | .00189601 | — |
| 1 | 14.09613549 | .00197335 | .00278858 | .25 |
| 2 | 8.08432634 | .00180812 | .00382762 | .5 |
| 3 | 7.17935248 | .00171677 | .00332277 | .25 |
| 4 | 6.19709301 | .00168036 | .00328348 | .25 |
| 5 | 6.17376652 | .00166284 | .00337256 | .25 |
| 6 | 6.04955643 | .00166231 | .00338964 | .0625 |
| 7 | 5.96706216 | .00166047 | .00342167 | .125 |
| 8 | 5.17313059 | .00163793 | .00353566 | .5 |
| 9 | 5.13728992 | .00163787 | .00354270 | .0625 |
| 10 | 5.07435571 | .00163763 | .00355586 | .125 |
| 11 | 5.03976638 | .00163762 | .00355887 | .03125 |
| 12 | 4.99549235 | .00163758 | .00356472 | .0625 |

The single run took 252.76 s and reached the inherited 12-update limit.
There were 45 trials, with 12 accepted and 33 rejected. All 13 completed
iterates have durable state, residual, constraints and history. Twelve
current matrices were published. Assembly, factorization, direction
calculation and line search took 114.57 s, 2.32 s, 28.98 s and 44.62 s
respectively. The residual dropped 68.9% below the R3 plateau with a
sign-resolved positive wall trace throughout.

The terminal coupled maximum is 4.99549 in the k block. The two momentum
blocks are .25338 and .22149, the epsilon block is .0035392 and the
constraint block is .0039919. The root criterion fails. The wall minimum
occurs at x=1.56414546 and the existing finite scalar inventory passes.
The covariance minima are .000437899 and .177077631. The fixed region
fails on the outer gradient departure, .00356472 at
`(1.6666666666666665,0,2)` in the outer-side `G[x,z]) component of
element 34. The minimum epsilon ratio, .974402, also falls below .98.
The k and sampled-K ratios remain within their intervals. The original
h0 and integrator history slots remain at initialization.

All 45 trial states passed the finite scalar and ground-trace checks.
The 33 rejections were residual non-descent, with a maximum of six trials
in any update. Lambda stayed one on the completed path. The pilot reached
its update limit without re-entering the fixed outer box.

The independent residual maximum is 4.995492353429094, and the maximum
scaled block of the native-to-independent difference is 2.96e-12.
Independent scalar reconstruction had zero unresolved signs. The terminal
ground polynomial was separately reconstructed at 110 digits and agrees
with the native minimum. The integrator's original history slots and
durable residual bytes were crosschecked. No native startup step was
authorized by the root and applicability results.

| Terminal independent observation | Value |
| --- | ---: |
| Weak continuity, inner / outer, scaled | 1.61e-15 / 1.67e-15 |
| Gauge, scaled | 2.17e-17 |
| Strong divergence L2, inner / outer | .0340962 / 8.31291e-5 |
| Normal / tangential total-traction jump RMS | 3.60649e-4 / 3.69812e-4 |
| k-flux jump RMS | .00130735 |
| Source accounting E minus D minus C | 4.50e-17 |
| Actual-state energy weak work minus direct work | 5.15e-14 |
| Free energy work at the terminal iterate | .0435202 |

The R4 energy diagnostic uses the external lift and the retained matching
epsilon coefficients for both terms. The native weak-minus-direct value is
-1.23e-12, and the native/independent energy crosscheck passes at 1.29e-12
against 3e-10. This is a same-state accounting check at the retained
non-root. The traction defects and the nonzero constraint block still
describe that iterate.

The complete local and pilot receipts live under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r4/epoch-003/cases-001`.
The pilot progress stream and its atomic last snapshot bind the complete
terminal payload. The saved-state outer audit and unsuccessful local
algorithm attempts remain in epochs 001 and 002. The repository holds the
compact [decision ledger](cls-epawf-revised-evidence/m8-r4/decision-ledger.json).
The full matrices, field vectors, wall extrema and independent audit are
retained on Linux.

## Next decision

At the retained pilot terminal state, test one coupled direction that
restores the existing outer-box inequalities alongside scalar feasibility.
The pilot cleared the observed wall-contact barrier but left both an outer
gradient excess and an epsilon-ratio deficit. A fixed-state comparison
can measure whether actual residual descent and recovery of the original
certificate can occur together before committing to another stage solve.

M8 remains open. Curved-wall MMS, active-cap curved and variable-H coupling,
fully discrete positivity, refinement, long-time robustness, scalability,
and physical floor/hill validation remain subsequent gates.
