# CV3-R1 numerical control: preconditioned-residual Newton map

This optional numerical change is evaluated before selecting a faster R1
development recipe. It changes neither a turbulence coefficient nor a wall
equation. Integrated strip continuity, the pressure multiplier, dynamic wall
epsilon and the independent 1e-7 physical residual/budget gates are retained.

## Fixed-point equivalence

For a transported block, let the current assembled free-row equation be

`A(x) y = b(x)` with nonlinear Dirichlet trace `y_D = g_D(x)`.

The existing Picard map solves this current system to an inner tolerance for
every Newton residual and finite-difference probe. Instead, project the trial
into the CG space, impose `g_D(x)`, and compute the **current** residual

`r_F = b_F(x) - A_FF(x) x_F - A_FD(x) g_D(x)`.

Apply the cached positive native inverse P once, giving

`F_D(x)=g_D(x)`, `F_F(x)=x_F + P_FF^{-1} r_F`.

For a nonsingular P, `F(x)=x` requires the same Dirichlet trace and `r_F=0`.
No current diffusion, mass, grad-div, Robin/normal coupling, source, convection,
nonlinear stress or boundary load is frozen in r. Only P retains its cached
coefficients. Projecting an off-space quadrature trial defines the numerical
extension used by Newton; physical fixed points are CG fields.

The existing pressure-defect and wall-reaction updates are retained. At a
fixed point their increments vanish, imposing the same weak divergence and
projected integrated-strip constraint. The nonzero outer relaxations do not
change these fixed points. Accepted fields retain strict positivity; a signed
linear map image is not a closure input or an accepted physical state.

`RANSNewtonMap=CurrentSolve` remains the default. The optional
`PreconditionedResidual` value is restricted to terrain Newton--Krylov with
the matrix-free operator. History and observation provenance bind the mode.
The reported matrix-free iteration count includes one Richardson sweep per
transported block in this mode; it does not mean a current Helmholtz system
was solved to convergence inside that sweep.

`RANSMatrixFreeTolerance` controls the inner current-operator solve on the
default `CurrentSolve` path. The residual map instead uses the native inverse
tolerance for its single sweep; acceptance still depends on the independently
assembled physical residuals. Its input retains the current-operator setting
for explicit provenance, rather than claiming an inner solve to that target.

## Validation contract, declared before new-mode runs

- Build in the private Release Linux tree, with an immutable candidate binary
  set; leave the running CurrentSolve binary intact.
- Reproduce accepted unforced small-step standard/CLS composite solutions,
  including S/R, using identical physical inputs, mesh and timestep. Require
  independently audited free/weak/projected/mass errors <=1e-7, positivity,
  exact prescribed correction coverage and independent strip inventories.
  Compare native complete fields to CurrentSolve within 1e-6 scaled error.
- Exercise representative full-history continuation/rank changes and mode
  mismatch rejection. A mode-specific history cannot silently change maps.
- Repeat registered tests and the current-solver reference/recovery controls
  before any R1 closure claim with a new numerical binary.
- Select a faster map only if its measured runtime improves and its physical
  residual/field checks pass. Failure remains a numerical control result;
  it does not establish model nonexistence. No failed stages enter statistics.

This is a solver-preconditioning control, not a closure regularization,
impermeable-wall replacement, temporal/spatial convergence result or
calibration-data selection.

## Pointwise scalar scaling control

The experimental epsilon field spans about five orders of magnitude. The
existing `FieldMaximum` Newton scaling weights an entire scalar by its
maximum; positivity at its smallest value can then restrict a physical
finite-difference step far below the requested field-relative step. This is
a conditioning hypothesis, not a demonstrated cause of the flow rejections.

Declare `RANSNewtonScalarScaling=Pointwise` as an optional numerical control,
retaining `FieldMaximum` as the default. Freeze positive reference values
`s_k(q)=k_incoming(q)` and `s_e(q)=epsilon_incoming(q)` for each physical
stage. Use `volume_metric(q)/s(q)^2` in the scalar Newton merit/Krylov norm,
and bound each finite-difference scalar direction relative to that same local
reference. Velocity, pressure and reaction scales are retained. Physical
unknowns, signed linear images, strict accepted-state positivity and all
independent equation/continuity gates are retained.

For positive finite s this is nonsingular row/column scaling of the same
Newton residual, with a different finite-difference step. Its zeros remain
the same assembled equations; no scalar is transformed, reset, limited or
frozen in those equations. Reference scales change only between physical
stages, and are not extra transported quantities. Bind the mode in history
and observation provenance, rejecting accidental mode changes.

Before selection, compile an immutable candidate, repeat the default Release
suite and all current-binary analytical controls, compare accepted S/R
standard/CLS solutions against the recorded field-maximum map within 1e-6,
and verify representative history/rank changes and scaling-mode mismatch
rejection. Test the refined-ramp first stage and its recomputed linear action
residual explicitly. Only fresh, fully developed four-case trajectories can
close R1; improved iteration counts or a short matching checkpoint cannot.
