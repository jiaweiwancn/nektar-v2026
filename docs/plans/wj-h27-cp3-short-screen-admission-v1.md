# H27 CP3 bounded cost and startup screens

Use the exact source19 application, table, Linux host and archived 700-cell
source19 startup session. This work can measure short startup viability and
the conventional compute comparison while the developed reference gate is
open. Positive recipe selection still needs the qualified reference and
the frozen finite perturbations, h/p/q and operating cadence.

The retained H26 conventional comparison is GLL8×8 with a verification
forcing, and H25R1 S0 used a different integration/output recipe. Reuse
their qualification context. For the cost pair, rerun one 12-step source19
startup at .005 and three serial 12-step stock VCSWeakPressure blocks on
the same physical mesh and executable, with tensor Gauss–Legendre 8×12
governing quadrature per stock cell. The v2 coupling uses that order per
geometric/material panel. Stock has no material panels. Keep Q3 velocity,
Q2 continuous pressure, IMEX(2), physical molecular viscosity, the same
historical velocity trace formula, and the same data/velocity boundary
types where the two formulations permit. Stock uses its own projection of
the analytic profile and high-order pressure Neumann boundaries except for
zero Dirichlet pressure at the outlet. It has zero componentwise natural
velocity data instead of v2 total traction, no model stress or wall loads,
and no transported scalars. Set NoCollection explicitly and start with SVV
disabled. Initial state projections and the pressure equations differ.

Both paths print each step and write every third step plus final fields.
Stock writes its own three-field snapshots. The coupled process additionally
exports complete restart and scientific data and rotates recovery writes.
Include all this work in the whole-process comparison and record file
counts/bytes. Stock's own internal step timer excludes checkpoints and some
post-step work, so keep it separate from complete wall time. Use [4–6],
[7–9] and [10–12] as warm step blocks. Accept a stock control only on
successful finite advancement to .06 and actual on-cadence outputs.

For the startup accuracy screen, use dt .02/.04/.08 to the common time .24
from the unchanged retained-inlet initialization. Write at .08/.16/.24.
Use the existing inverse policy, maximum age eight, relative drift .02 and
tolerances 3e-12/1e-12. At .04 separately test age24 and nonlinear1e-10,
each with every other numerical input fixed. Compare matched accepted
fields with the physical owner-aware norms at the common endpoint. Retain
the current independent momentum, continuity, k and epsilon checks and
original physical scales. A failure in one temporal/controls path does
not stop independent bounded cases. The startup clock is shorter than the
required pilot and has no developed perturbation. Results are a screen.

Charge fixtures 30 s, one paired startup 130 s, three stock blocks 45 s
each, the temporal runs 180/120/90 s and two .04 controls 120 s each. Admit
four independent current-equation/common-time assessments at up to 160 s
each, a 50 s cost assessment, and 90 s publication. Total cap 1705 s
under the live 43,200 aggregate ceiling with 3,600 protected for evidence.
Check the complete cap, host allocation and space after CP2 operations
settle. No build, one rank, six thread limits one, 72 actual workers, science
RSS and both storage caps 8 GiB and free reserve 15 GiB. Use create-only
paths and preserve failed receipts. No long physical pilot or full Task 3
horizon is permitted by this screen.

Admit a separate create-only assessment for up to 50 seconds. The first
summarizer attempted to read the six-component residual vector as a mapping
and exited before publishing. It also selected the Picard per-block linear
counter. Re-read the exact native step logs using the residual vector and
active coupled Newton–Krylov counters. Report conservative exclusive
whole-process attribution. Preserve the failed receipt and all timing and
scientific evidence. No native run needs repetition. The adjusted phase
cap is 1,755 seconds.
