# WJ-H27 native steady preparation

The opt-in `WJNativeSolve=SteadyReference` mode solves the current native flat
two-dimensional equations at zero physical time. It requires `NumSteps=0` and
an explicit `WJSteadySeconds` ceiling. It does not replace the retained inlet,
change the physical boundary laws, or interpret nonlinear iterations as
physical spinup. A converged result receives fresh step-zero current/previous
states and freshly evaluated explicit histories. A failed solve preserves the
input history and writes a `NOT_CONVERGED` diagnostic.

The unknown vector contains each actual component's free velocity coefficients,
mechanical pressure, the mean constraint where its actual nullspace requires
one, free full-domain k, and free outer epsilon. At every residual evaluation,
the epsilon interface coefficients are reconstructed from `A(K)/H`. Current
advection, full stress, dependent dissipation/cutoff, scalar diffusion and
sources are refreshed. The equation residual contains no time-mass term.

The method is feasible Newton with flexible right-preconditioned GMRES and
second-order finite-difference Jacobian actions. Both signs of each trial use
the complete current equations. Only an explicitly typed inadmissible state
may select a feasible one-sided stencil or reject a line-search trial; input
and programming errors remain fatal. The norm is the fixed physical mass-dual
block maximum already used for the coupled equations. There is no large-alpha
normalization. The final acceptance also checks physical continuity after
eliminating the gauge multiplier.

The inverse approximation is rebuilt at each Newton base. Its flow block uses
current symmetric-strain diffusion, the actual pressure/AWF action and the
current momentum-convection derivative. Scalar blocks use current diffusion,
conservative transport and positive epsilon destruction. The triangular action
includes the actual Q3 known-k force and the derivative of the projected
dependent epsilon trace. An explicitly declared nonnegative
`WJSteadyPreconditionerMass` may improve this inverse approximation; it never
enters the actual steady residual or its Jacobian action.

For the verification-only frozen MMS, analytic u_t, w_t, k_t and epsilon_t
are subtracted from the existing continuous manufactured forcing at t=0.
The exact frozen fields therefore solve the continuous steady equations.
The independent oracle must make the same continuous-time specialization
independently; no projected/discrete equilibrium source is constructed.

Qualification requires the standalone analytic nonlinear/Jacobian/error-control
probe, an independently checked native coupled steady MMS, and the actual
same-inlet physical reference plus a finer calculation. None of these are
inferred from another. The next bounded stage admits an immutable four-worker
build up to 600 seconds, focused helper operations up to 90 seconds each,
independent tiny checks up to 300 seconds, and a physical one-step diagnostic
up to 180 seconds. A physical steady diagnostic may be admitted only after the
native steady-MMS checks pass, with an explicit prefix cap and a new forecast
before a larger reference/refinement matrix.

Exact immutable native H traces are evaluated from the validated Modified_A
edge polynomial and its actual signed scalar/component maps. Construction
compares value and physical tangential gradient against native owned basis
evaluation on nonconstant synthetic coefficient fields and five points per
H edge, including both endpoints. No current field value is retained across
state changes, and distinct native gradient owners are not merged.
