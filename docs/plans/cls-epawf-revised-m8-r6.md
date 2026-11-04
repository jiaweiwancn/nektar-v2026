# M8-R6 coupled restoration direction accuracy

M8-R6 uses the retained R5 iterate 1 with its own current J1. The original
initialized h0 remains the history for alpha=400 at time .0025, eta=0 and
zero prescribed forcing. SaturatedR2 with rho=.8, beta=.4, GammaScale=1,
GammaUniform=0, the original outer response, viscosity and scalar diffusion
conventions, the projected epsilon trace and capped source with theta=.2
and delta=.1 remain as in R5. The outer6 p3/p1 spaces and quadrature give
2,165 unknowns under the 3,000-unknown limit. The coupled root gate is
2e-12 and the independent action gate is 1e-7.

The reviewed revision was `7ae3daae99c1f3e8fb9be1ea357ef01ef9169293`.
The baseline replay loaded J1, verified its association with the saved
state, residual and h0, and reproduced the retained constraint bytes.
The state has (E=3.92879913247) and (V=1.92353785366). The retained
gradient departure is .00292353785 and the minimum retained epsilon ratio
is .97580135906. Those are actual state observations. The cut violations
below are defects in the search direction's linear model.

## Reproduced numerical failure

| Final entering inequality | Full restoration | Partial restoration |
| --- | ---: | ---: |
| Native gap | -3.57250e-11 | -9.01751e-12 |
| Independent coefficient-space gap | -3.57251e-11 | -9.01750e-12 |
| Independent residual-space gap | -4.32e-16 | +1.96e-16 |
| Transpose-solve contribution to gap difference | -3.16970e-11 | -1.01429e-12 |
| Forward reconstruction contribution | -4.02401e-12 | -8.00377e-12 |
| Consistent-mass transform contribution | -3.82e-15 | +3.61e-16 |
| Exact recurrent entries across 96 scans | 54 | 71 |
| Active normalized-row condition number | 11.10 | 27.30 |

Both last entering cuts are velocity (x,z) gradient inequalities.
The full restoration witness is at (-3, 0, 2.3333333333) in element 30.
The partial witness is at (-1.5833333333, 0, 2) in element 32 on the outer
side of the matching interface. The native and independently accumulated
coefficient gaps agree closely. Each entering row has an exactly equal
coefficient vector and right-hand side already in the active set. The
replay released and re-added this row while the physical gap remained.
The final active sets had 29 and 14 rows and minimum QR diagonals .2171
and .07851. The corresponding multipliers were dual feasible,
stationarity was below 1e-15, and complementarity below 7.5e-16.

For a selected row, let (x) be the computed transpose solve, (w=d-d_N),
and (n) the unscaled computed normal. The measured gap difference
decomposes into
[
  (a-J^T x)^T w + x^T(Jw-Bz) + (x^TBz-n^Tz)
]
plus rounding of the affine shift. The two solve terms explain the
observed disagreement. The full and partial transpose defects in
infinity norm were 6.84e-10 and 1.66e-10, and the forward correction
defects were 2.32e-15 and 1.67e-15 before contraction with (x). The
complete (Jd+R-Bz) defects were 1.19e-15 and 2.89e-15 in infinity norm.
The retained algebraic residuals were below 3.2e-16 in coefficient units.
The 60-digit active problem with the saved transformed normals changes
the largest (z) coefficient by less than 3.4e-16 and keeps the
re-entering cut at transformed equality. It checks the small active
solve. The separate gap decomposition checks the transforms that
produced its input.

## Repair and trial policy

`RANSWallDirectionAccuracy=RefinedR6` opts in alongside
`JointRestorationR5`. `ReplayR6` is the diagnostic-only legacy
arithmetic replay on a bound matrix. R6 retains the original pivoted LU,
uses two extended-residual iterative refinement corrections for the
forward and transpose solves, and accumulates the consistent-mass
products in extended precision. The recovered direction solves directly
against (Bz-R). The active set still uses the same residual-space QR,
physical consistent mass, full algebraic restoration, dimensionless
radius and scalar and ground safeguards.

The full original inventory is scanned in native and separate extended
row arithmetic. Both scans must meet the original 5e-12 direction
criterion. If an already-active identical row still violates it, the
solver records the failed candidate and classifies the direction
accuracy as unresolved. No inequality is removed, and there is no
inexact acceptance or raised tolerance. Actual trials still rebuild
the scalar inventory, ground extrema, dependent trace and outer box.
The residual and outer acceptance rules retain the R5 thresholds and
the per-family nonincrease rule.

A continuation update reserves 16 trials for full restoration and 16
for partial restoration. The bounded fixed-state search uses four for
each. Partial restoration remains available outside the box after
either a full direction failure or a completed full line search with
no acceptable trial. The partial model still targets a 25% decrease in
each positive family maximum. There is no change to the final
qualification limits.

## Fixed-state verification

| Check | Saved R5 terminal state | Earlier R5 partial-restoration control state |
| --- | ---: | ---: |
| Full-direction scans | 50 | 51 |
| Worst native cut violation | 6.514e-13 | 3.399e-12 |
| Worst extended original-row violation | 6.514e-13 | 3.399e-12 |
| Largest native-to-extended scan difference | 2.42e-16 | 4.38e-16 |
| Largest retained coefficient-to-transformed gap difference | 2.98e-15 | 2.82e-15 |
| Complete (Jd+R-Bz) coefficient infinity defect | 2.25e-18 | 2.35e-17 |

The refined full solve completed with the original criterion in both
cases. The control is the state from which R5 needed its partial
restoration. R6 constructed a full direction there and reduced the
native residual to 3.98843232, with gradient excess 4.997e-12 and no
epsilon excess on the actual trial. The gradient result remains strictly
outside the fixed box.

At the saved terminal state, the refined full trial reduced the maximum
outer excess to 1.280e-12 but increased (E) to 7.91068283. The half
trial passed the existing joint acceptance gate. It reduced (E) from
3.92879913 to 3.72003500 and the gradient and epsilon excesses from
1.92353785 and .20993205 to .96176893 and .10496602. The ground minimum
was .0016286283. The scalar inventory passed, lambda stayed one and the
profile branch pattern was unchanged. The actual constraint block was
.00136244 at the damped trial, separate from the direction's linearized
algebraic restoration.

The independent R3 value oracle checked the actual retained full
direction with the same bound J1 at complex steps 1e-20 and 1e-30.
The maximum physical action error was 1.544e-14. The largest normalized
derivative error on the governing outer-box witnesses was 6.51e-16.
The checks treat the retained epsilon polynomial and projected trace as
separate observations.

The refined direction at the terminal state took 17.81 s. The
continuation budget includes 20% headroom on the prior measured current-J
assembly, the refined direction cost, six trials at the local mean
cost, publication and a 90 s setup and terminal audit allowance. The
estimated cost was 46.513 s per update, with 17 permitted updates and a
900 s native cap. The full and partial searches share the 32-trial cap
through their separate allocations.

## Single exact-target continuation

The first update reproduced the fixed-state half trial. The continuation
rebuilt its own current matrix at each successive accepted state. The
first rebuilt matrix matches the bound R5 J1 bytes. Seventeen direction
solves passed both complete original-row scans. The maximum native
and extended cut violations over the run were 3.272e-12, below 5e-12,
and the maximum number of cut scans was 54. The final published current
matrix belongs to iterate 16. Iterate 17 is the completed state saved
at the update limit.

| Completed update | (E) | (V) | Accepted fraction |
| ---: | ---: | ---: | ---: |
| 0 | 3.92879913 | 1.92353785 | — |
| 1 | 3.72003500 | .96176893 | .5 |
| 5 | 3.48853329 | .35666197 | .5 |
| 10 | 3.35484521 | .11520112 | .25 |
| 15 | 3.24588037 | .03100811 | .5 |
| 16 | 3.24091990 | .02713210 | .125 |
| 17 | 3.23071334 | .02374059 | .125 |

The native process used 585.27 s and reached the declared update limit
while both metrics were decreasing. There were 61 trials. Seventeen were
accepted and 44 increased the complete residual. The largest number of
trials in any update was six. All accepted updates used the full
restoration direction, damped when needed. At the last update the full
linear model predicted (E=3.12702) at fraction one. Its actual residual
was 11.66103. Fraction .125 gave actual (E=3.23071). The continuation
therefore tested rejection and damping of constructed directions. Its
full line searches never exhausted their reserved allowance.

| Separate terminal outcome | Result |
| --- | --- |
| Coupled root | Unresolved. The k block controls (E=3.23071334) against 2e-12 |
| Scalar evaluation domain | Passed on the original finite inventory |
| Continuous ground trace | Minimum k .00162860512, positive |
| Fixed constitutive region | Outside. Gradient departure .00102374059 and minimum epsilon ratio .97994817835 |
| Measured covariance | Minima .00045732251 inner and .17687429415 outer |
| Lambda | Minimum one |
| Physical history | Original h0 and startup slots retained, zero physical steps |

The terminal gradient witness is (G_{xz}) at
(-2.66666666667, 0, 2) on the outer side of element 31. The positive
gradient and epsilon family excesses are .02374059 and .00259108.
At (1.75, 0, 2) on the outer side of element 34, the retained epsilon
ratio is .97994817835 and the projected target ratio is .97994820597.
Both remain below .98. The two momentum residual blocks are .14793540 and .23717505, the
epsilon block is 4.74e-5 and the constraint block is .00142453.
The native weak-continuity blocks are 1.04e-15 and 8.62e-16 and the
gauge block is 3.41e-17. Incremental ellipticity remains unqualified at this out-of-box state.

The independently reconstructed residual maximum is 3.23071334389360.
The largest scaled native-to-independent block difference is 2.84e-12.
The independent checks corroborate the finite scalar domain, outer-box
status and measured covariance. The extended ground reconstruction
agrees with the native minimum. The actual epsilon constraint remains
nonzero. The source defect E minus D minus C is -1.72e-16 and the
same-state energy weak-minus-direct work is -3.82e-14. The native and
independent energy diagnostics agree to 1.06e-12. Independent normal
and tangential total-traction jump RMS values are 3.86e-4 and 2.05e-4.
The k-flux jump RMS is 6.54e-4. Strong divergence L2 is .03677 inner
and 8.48e-5 outer. The source and energy identities concern the retained
non-root, whose actual constraint and traction defects are reported
separately.

The diagnostic replay was built at `984abbd1db3b9bf0a6a8c9f57943a49c66e3f124`.
The refined experiment and continuation were built at
`7f637f98367a4e99cd0367d555f28452b3d0a27a`. Builds used 16 compiler
workers and the installed MPI Release libraries. The native solver used
one rank and OpenMP, OpenBLAS, MKL, BLIS, vecLib and NumExpr thread counts
of one. The terminal replay and continuation ran on CPU 0 and the
earlier-state control on CPU 48. The independent qualification and
terminal audit ran on CPU 49.

The three focused native time, evolution and wall-stage tests and nine
predictor binding tests passed. The fallback policy regression compiles
both triggers, exhaustion after a constructed full direction and an
unresolved full direction, as well as the 16+16 and 4+4 limits. That
revision passed the Linux wall-stage smoke test. The continuation keeps
its original producing source and executable. The failed first
diagnostic build and its compiler receipt remain in epoch 001.

The [decision ledger](cls-epawf-revised-evidence/m8-r6/decision-ledger.json)
pins the frozen matrices, stage inputs, current matrices and all terminal
results. The Linux evidence root is
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r6`.
Epoch 002 holds the legacy cut replay, epoch 003 the refined fixed-state
controls, actual trial checks, continuation and independent audit. The
complete row vectors, active histories and checkpoint payloads remain
outside the source checkout.

## Next decision

At the saved R6 terminal state, test one model-based trust-radius
contraction while retaining the same scalar, trace and outer-box
constraints. The cut solves now meet the original criterion, and the
last full trial's residual was 11.66103 against a model prediction of
3.12702. Smaller trial steps preserved joint descent. A local comparison
can establish whether a shorter constructed direction improves the
native residual reduction before budgeting any further continuation.

M8 remains open. Curved-wall MMS, active-cap curved or variable-height
coupling, full-element and fully discrete positivity, refinement,
long-time robustness, scalability and physical floor/hill validation
remain later gates.
