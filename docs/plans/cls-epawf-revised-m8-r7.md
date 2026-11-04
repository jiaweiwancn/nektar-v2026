# M8-R7 finite-step progress at the R6 endpoint

The reviewed source was `d0ac41af5b249d61ccc5a97a1f0db016a9567f19`.
The nonlinear guess is the retained R6 iterate 17. The exact target keeps
alpha=400, time .0025, the original h0 and startup slots, eta=0 and zero
prescribed forcing. SaturatedR2 with rho=.8, beta=.4, GammaScale=1,
GammaUniform=0, the original outer response, viscosity and scalar
diffusivities, the projected epsilon trace and the capped source with
theta=.2 and delta=.1 are unchanged. The outer6 p3/p1 spaces and quadrature
give 2,165 unknowns. The root, independent action and linear cut
thresholds remain 2e-12, 1e-7 and 5e-12, respectively.

## Matched comparison protocol

`RANSWallR7Policy=Compare` opts into the fixed-state experiment with
`RefinedR6`. The two forward/transpose refinement corrections, direct
reconstruction from $Bz-R$, full linearized algebraic restoration and
both full original-row scans apply to every candidate. The native
comparison assembles its own J17. Its input contains the earlier J16 pin
for provenance and does not supply that matrix to the solver.

The last accepted R6 displacement has dimensionless infinity norm .25.
Its radius-2 direction had norm 2 and the accepted fraction was .125.
The first contracted radius is the smaller of half the measured J17
radius-2 direction norm and twice the last accepted displacement. This
puts the re-solved direction near the previously observed finite-step
regime while giving it room to choose a different direction. One
adjustment is reserved for a failed construction or a nonbinding radius.

The radius-2 control has eight local trial slots. Each contracted radius
has four for full restoration, and four for the inherited 25% partial
outer restoration if full construction fails. The single guard-only
control, when indicated, has eight slots. It uses the same radius and
outer model as its guarded counterpart. It removes nondepletion only at
strictly positive ground values and retains the 5% reserve, contact
handling and continuous ground extrema. A summed multiplier-times-base-k
quantity greater than 1% of the quadratic model objective triggers this
test. This is a screening measure in the existing normalized cut units.
The actual trial comparison decides whether the guard change is useful.

All local candidates use the continuation acceptance requirement
$10^{-4}f$ at trial fraction $f$, including per-family box
nonincrease. The actual nonlinear trace residual is recomputed in the
normal coupled residual. Predicted
and actual reduction are measured in the same maximum scaled residual
norm $E$. The agreement ratio is recorded when predicted reduction is
positive. Exhausted trials and unresolved direction constructions have
distinct records.

## J17 comparison

The terminal residual, state, history and constraints reproduce the R6
iterate 17 payloads byte for byte. The base has $E=3.23071334389$,
$V=.02374058719$, ground minimum .00162860512 and a passing original
finite scalar inventory. A fresh J17 took 9.197 s to assemble and .193 s
to factor. Its hash differs from J16. All candidates below used this
same current matrix.

| Model | Radius | Direction norm | Full predicted $E$ | Full actual $E$ | Accepted fraction | Accepted $E$ | Accepted $V$ | Agreement at acceptance |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| R6 control | 2 | 1.07552 | 3.12701 | 3.54371 | .25 | 3.22279 | .0178054 | .17618 |
| Contracted, guarded | .5 | .5 | 3.12701 | 3.71053 | .25 | 3.22457 | .0178054 | .13649 |
| Contracted, positive ground control | .5 | .5 | 2.91761 | 3.24589 | .5 | 3.12896 | .0118703 | .56109 |

The first contraction binds and the full algebraic and outer
restorations are both constructible. No second radius was needed. The
necessary algebraic radius estimate in the ledger is a rowwise bound
on the dimensionless coefficient infinity norm. It is not a sufficient
feasibility test. The smaller radius alone produced less actual
reduction with the same number of trials. The full-step model error
also grew. In the guarded solves the positive-ground screen was .1413
of the model objective. The matched guard-only control improved the
full model, accepted progress and agreement. Its accepted ground
minimum was .000855645 and its rejected full trial still retained a
positive ground minimum of .0000815492. The positive-ground control
provides the evidence for changing this numerical restriction at the
tested radius. The radius-2 result supplies the ordinary damping control.

Direction and line-search times in seconds were 13.54 and 5.79 for the
control, 13.18 and 5.90 for the guarded contraction, and 17.75 and 3.88
for the positive-ground control.

The ground, scalar, covariance and lambda observations, actual and
predicted residual block norms and outer-family violations, trace
residuals and active row slacks are retained separately for every
evaluable trial. Full accepted-direction algebraic restoration is
checked independently of the actual damped trace residual.

The positive-ground row is at (1.56236145442, 0, 0). Its multiplier
in original normalized cut units is 218.110 for radius 2 and 218.107
for radius .5. The shorter direction also activates the bound on
coefficient 1536. The guard-only direction activates the continuous
5% ground reserve at (1.56365327140, 0, 0) with multiplier 202.603.
The native slacks of these active rows are within 1.4e-14 of zero.
The necessary algebraic radius lower bound is .000135700. The three
directions' full-model algebraic blocks are below 9.4e-19.

## Finite box entry and selected numerical policy

The accepted boundary-targeting guard-only half step left about half
the excess in each failing family. For a linearized signed bound with
base excess at most $e$ and full-step target $-m$, fraction $f$
has excess at most $(1-f)e-fm$. The matched interior probe reuses J17
and keeps the .5 radius and positive-ground control. It fixes
$m=1.25e$ in the initially failing gradient and epsilon families,
with a numerical minimum of 1e-8 in such families. The actual margins
are .02967573398 and .00323885299 in the normalized cut units. At
$f=.5$, the affine bound has clearance $.125e$, far above the
5e-12 linear cut threshold. The actual projected epsilon and the
retained epsilon are tested separately. The full trial entered the box
and increased $E$ to 3.24582. The accepted half trial also entered the
box and lowered $E$ to 3.12912276, with agreement .56086.

The pilot starts at the tested .5 radius with the positive-value
nondepletion restriction removed and the fixed search-only interior
target until actual box entry. The unchanged box test then maintains
membership. Positive excess, including arithmetic-scale excess, is
outside. The true-contact rule, the scalar and continuous-wall reserves
and all original nonlinear scalar and physical tests remain in effect.

The radius grows by a factor of 1.5 up to 2 after a full step with
agreement at least .75. It halves down to .125 after agreement below
.25 or acceptance at a fraction at most .25. Otherwise it stays fixed.
If the first full search fails, one further full construction is
available. It contracts after nonlinear trial failure or expands
after construction failure, within the stated radius limits. Partial
outer restoration keeps its 25% target and full algebraic restoration.
The full searches share 16 trials as 8+8, and partial restoration
reserves 16. There are at most three constructions at each base state.
Every accepted state triggers a new current Jacobian. The hard limits
are 32 accepted updates and 900 native seconds, with no restart.

The selected J17 action passes the independent 1e-7 test with maximum
physical action error 4.03e-15. The largest normalized outer-witness
action difference is 1.64e-15. The check uses the same J17 bytes and
the unchanged value oracle. The interior direction, factorization and
trials took 20.26 s, in addition to the measured J17 assembly. The
continuation reserves the hard limits for changing per-state cost
instead of using R6's 17-update cap.

The comparison was produced at
`bc607ec3ba0cf852028d9516d968d586abc4776a` and the interior probe and
continuation at `a5e3aa2b4a060fda02be113e144f3bc86e644338`.
Both isolated Linux Release builds used 16 compiler workers. Native
stage runs used one MPI rank on CPU 0 with OpenMP, OpenBLAS, MKL, BLIS,
vecLib and NumExpr thread counts of one. Independent checks used CPU 49.
The shared Linux checkout and installation were preserved. The three
focused native time, evolution and wall-stage tests passed on both
epochs. The nine predictor binding tests passed on the continuation
epoch.

## Single exact-target continuation

The first current matrix matches J17 byte for byte. The first update
reproduces the accepted interior-probe half step and enters the original
box. Every later accepted state remains inside. The continuation reaches
the 900 s native wall cap after 25 accepted updates. The supervisor
records exit 124 and 900.18 s including termination. The 32-update
allowance is not exhausted. The terminal iterate has its own complete
current matrix. A direction was built from that matrix before the cap.

| Accepted update | $E$ | $V$ | Ground minimum | Accepted fraction |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 3.23071334 | .02374059 | .001628605 | — |
| 1 | 3.12912276 | 0 | .000855647 | .5 |
| 5 | 3.06655755 | 0 | .000446822 | .125 |
| 10 | 2.96164006 | 0 | .000105976 | .25 |
| 15 | 2.93817480 | 0 | 4.08085e-5 | .125 |
| 20 | 2.93316099 | 0 | 2.63177e-5 | .0625 |
| 24 | 2.92803060 | 0 | 1.88979e-5 | .03125 |
| 25 | 2.92612194 | 0 | 1.44097e-5 | .25 |

Of 110 recorded trials, 25 were accepted, 84 failed the actual residual
decrease rule and one failed the actual outer-family rule. There were
at most seven trials and one direction construction per update. The
full retry and partial allocations remained available. No accepted
state failed the scalar or continuous ground tests. Rechecking the
durable records confirms the acceptance, radius-update and work limits.
The largest original-row cut defect across both scans was 4.954e-12,
within 5e-12.

At the terminal state the k block dominates $E$. The momentum blocks
are .11851117 and .22035743, the epsilon block is 5.12e-5 and the
constraint block is .000867741. The weak-continuity and gauge blocks
are of order 1e-15 or smaller. The current direction at radius .125
predicts full-step $E=2.90370$. The continuous-wall reserve remains
active at (1.56385178842, 0, 0), with multiplier 206.516 in original
normalized cut units. A coefficient radius bound also remains active,
with multiplier 2.875e-5 in its normalized units. These are active
constraints of the local search model. The observed stop is the wall
cap during continued residual reduction.

The terminal direction's full trial gave $E=10.04468$ against the model
prediction of 2.90370. The .5, .25 and .125 trials also increased $E$.
The wall cap interrupted the search before the remaining local
allocation was exhausted.

## Independent terminal qualification

| Separate outcome | Terminal observation |
| --- | --- |
| Nonlinear root | The independent maximum residual is 2.92612194240211 against 2e-12 |
| Actual box and finite constitutive qualification | Passed on the original inventory. Gradient departure .0009999615834 and minimum epsilon ratio .980000083857 |
| Finite scalar positivity | Passed with zero unresolved signs in the independent inventory |
| Continuous ground trace | Positive minimum 1.4409725e-5, confirmed by independent 80- and 110-digit polynomial reconstruction |
| Covariance and lambda | Measured minima .00012084919 inner and .17687695302 outer. Minimum lambda one |
| Direction accuracy | Both original-row scans stayed within 5e-12. The representative J17 action passed 1e-7 |
| Nonlinear model agreement | The last full model substantially underpredicted the actual residual |
| Stop and physical history | Wall cap during the terminal line search. Zero accepted physical steps |

The native-to-independent maximum scaled block difference is 3.00e-12.
The actual matching epsilon trace differs from its independently
projected target by a maximum coefficient error of 6.30e-5, or
.000867741 relative to the epsilon scale, across 16 matching rows.
The reconstructed constraint targets agree with the saved targets to
4.01e-16. The source defect E minus D minus C is -1.43e-16. The
same-state weak-minus-direct energy work is -4.44e-14. The normal and
tangential total-traction jump RMS values are 3.83e-4 and 2.06e-4.
The k-flux jump RMS is 6.53e-4. The accounting identities concern
the saved non-root.

The supervisor stopped the native process before it published a final
native report. The terminal audit reconstructs the saved state
independently and uses the durable residual, constraint and history
payloads. It therefore has independent source and energy observations
and the original history binding, with final native budget tables
unavailable for crosscheck. The same-executable startup probe verifies
the original four history slots and zero startup steps. Every saved
continuation state retains the same h0 bytes with history advancement
false.

The [decision ledger](cls-epawf-revised-evidence/m8-r7/decision-ledger.json)
pins the local and continuation runtimes, snapshots, J17 and current
terminal matrix, trial and action reports, and independent audit. The
full numerical payloads remain under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r7`.

## Next decision

At the saved in-box iterate 25, compare one curvature-aware trial
predictor for the dominant k residual against ordinary damping of the
retained J25 direction. Keep the radius, continuous-wall reserve and
actual feasibility gates fixed. The terminal search gives a direct
test case because the linear model predicted reduction while all four
completed trial fractions increased the residual.

M8 remains open. Curved-wall MMS, active-cap curved or variable-height
coupling, full-element and fully discrete positivity, refinement,
long-time robustness, scalability and physical floor/hill validation
remain outstanding.
