# H27 current nonlinear iteration and fresh reference fields

The physical 700-cell diagnostic at source-v8 completed no physical step:
after 60 sweeps its current momentum residual was 1.585304203902587e-9,
above the unchanged 1e-10 native coupled target. Its k and epsilon residuals
were 2.0712428840603278e-12 and 5.7530637561806853e-12. Positivity and covariance
were checked during these sweeps; the final principal check was not reached.
This is nonlinear-solver evidence, not a new physical reference or cost pass.

The next opt-in numerical method is Anderson acceleration with depth 4
(default remains 0). It stores the current step's fixed-point images F and
scaled updates g=diag(s)^-1(F(x)-x). Differences form DG and DF; pivoted QR
solves min ||g-DG gamma||. The proposal is F-DF gamma. Scales are Uref=1,
pressure tau0, k Kin, and the fixed local epsilon reference RMS from its
weighted-mass diagonal. A relative QR rank threshold 1e-10 and coefficient
norm ceiling 1000 bound numerical extrapolation. They change no physical field
or equation. Fractions 1, 1/2, 1/4, 1/8 are tested, with fresh dependent Gamma,
cutoff, constitutive/scalar data and full coupled residual. A proposal is used
only if it passes the interim scientific guards and reduces that true residual.
Final acceptance still requires all original residual and principal margins.

Accepted histories, physical dt and EX2 states are untouched by provisional
iterations. Only a typed recoverable linear failure permits the one same-dt
cold-factor retry; exhausting nonlinear sweeps no longer repeats an identical
attempt from the old state. Diagnostic traces are opt-in. Internal failures
persist the last accepted complete history, the failed-attempt residual trace,
and their distinct provenance. A bounded steady failure additionally retains
its best evaluated unknowns/raw residual, finite-difference and linear-check
records, and a separately labelled admissible candidate for continuation.

`WJSteadySeed` is verification-only and restricted to manufactured steady mode.
It perturbs free scalar coefficients by the declared small reference-scale
amplitude and reconstructs a fully admissible fresh history. This exercises
the actual coupled Newton/Krylov/trace-elimination path; the previous polynomial
steady MMS was already exact at its initial state and exercised no Newton step.

`WJInitialFields` is separate from `WJRestart`. It imports step-zero, pinned
source/model/reference/geometry coefficient packets. Identical inventories
copy exact coefficients; a different supported native grid projects values
with explicit original owners and the actual signed maps. Current Gamma,
advection and both fresh histories are rebuilt using the new dt/rules. A
candidate is allowed only for steady preparation. Advancement from reference
fields also checks current stationarity on its selected recipe. The packet's
native steady-convergence role is not the independent physical/refinement gate.
Cold restart continues to require the complete existing history signature.
An explicitly pinned different `WJInitialSourceIdentity` is accepted solely as
a steady-preparation guess, with the same model/reference/geometry bindings.
The current equations, admissibility and full steady residual are re-solved
before it can become a current reference. It cannot bypass advancement or
cold-restart source checks.

The prospective next stage admits one immutable four-worker build up to
600 seconds, initial-field controls up to 90 seconds, coupled tiny MMS and
independent checks up to 330 seconds per operation, a 180-second accelerated
physical one-step diagnostic, and a 240-second unaccelerated control. A physical
steady diagnostic is conditional on the seeded native and independent MMS
checks, with at most 240 seconds inside Krylov and a 360-second external cap
including construction, initialization, forcing and final checks. No long
physical pilot/refinement matrix or Task3 horizon is admitted by this note.
All operations share the original H27 account and protected evidence reserve.
