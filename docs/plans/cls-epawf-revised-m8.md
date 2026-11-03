# M8 — Robust variant qualification and conditional terrain validation

## Material Passport

- Mode: experiment-agent implementation, run and reproducibility validation.
- Authority: user approved the three-variant proposal and requested M8 execution.
- Start: 2026-11-03, baseline `22d5c930380a202e29bb2c58613e3763bc2aa372`.
- Branch: `codex/cls-epawf-revised`.
- Status: **IN PROGRESS; no robust or terrain qualification is claimed.**
- Numerical environment: isolated Linux run root through Tailscale
  `jay@100.105.247.110`; compiler workers 16, numerical library threads one.
  All concurrent compiler, MPI and test workers must stay below 72.

## Approved scope and model decision

The approval selects a research sequence: common cubic scaling plus bounded
inner destruction first, then a reference-preserving positive-increment
regularization. A changed invariant response remains an alternative if the
second candidate requires excessive diffusion. Original CLS and standard
k-epsilon remain controls. The frozen revised note, M0 contract and original
kernel formulas are retained; M8 has its own explicit model selection.

`RANSWallVariant=Original` retains the old outer kernel. `CubicScaled` applies
one prescribed beta to all active cubic terms; beta=0.4 is the primary
reference-local candidate, 0.45 its sensitivity case. `Regularized` adds
`-gamma(x)*(S-Sref)` to the Reynolds stress of that same scaled law. One stress
must supply momentum, production in both turbulence equations, and interface
traction. There is no extra addition of gamma to scalar diffusivities.

The prescribed viscosity is
`gamma(x)=GammaUniform+GammaScale*kappa*u_star*r`. Its coefficients do not
depend on the evolving state. `Sref` is the symmetric incoming flat shear,
`Gref=flow_incoming tensor normal_incoming*u_star/(kappa*r)`. The incoming
frame, physical normal distance convention and extension over terrain are
recorded. This is an explicit reference-dependent constitutive model, not
physical molecular viscosity or a changed nonlinear iteration.

The selected inner source uses the capped-normalized destruction construction
from the earlier note's remedies supplement. Initial declared parameters are
theta=0.2, delta=0.1, piecewise-linear limiter; a C2 quintic is an explicit
separate option. For pre-repair E,D with equal metric ray budgets B,
`w=min(1,k/(theta*K))`, `b=integral(J*D*w)/B`, `m=max(b,delta)`, and
`lambda=b/m`. Replace destruction by `Dpos=D*w/m`, dependent dissipation by
`Epos=lambda*E`, and total redistribution by `Cpos=Epos-Dpos`. The net source
is actual production minus Dpos, assembled once. This same repair applies
to the selected chi=0 baseline. In severe depletion it deliberately reduces
the integrated loss and changes the dissipation trace. Negative represented
k remains invalid; no field floor or post-step clipping is introduced.

## Gates fixed before the campaign

1. Independent compiled-kernel identities, full indexed stress and branch
   derivatives; exact original identity; common-beta reference equality;
   production/energy consistency. Reproduce the note's reference interval
   certificate and the off-reference adverse shear witness independently.
2. Independent ray source quadrature and full directional dependence,
   plateau identity, zero-k reaction, rate bounds, capped/depleted regimes,
   actual trace and conserved redistribution. A continuum source identity
   does not imply a positive high-order CG/IMEX discretization.
3. A rigorous finite state-envelope tangent bound for the regularization,
   covering both limiting invariant branches. A sampled scan cannot stand
   in for this certificate. Record its required gamma and diffusion cost.
4. Native flat reference/stationary stages, full energy and coupled residuals,
   curved manufactured fields, interfaces, independent dt/p/mesh refinements,
   several perturbation amplitudes/wavelengths, MPI and complete restart.
   Match the existing scaled gates (residual 1e-9, isolated identities 1e-10,
   full-field MPI/restart 1e-8); use tighter internal nonlinear tolerances.
5. Separate covariance, outer k/epsilon and bounded source checks. Inner
   depletion controls test the source mechanism; they are not evidence of
   full coupled realizability. Long unforced trajectories must actually
   complete, with durations in advection and turnover units. Stopped prefixes
   cannot qualify transition passage or supply developed-flow averages.
6. Only after the preceding qualification: paired physical empty-floor/hill
   cases with the same repaired patch, reference inputs and baseline, with
   eta=0. Version the Cao-Tamura selector/freezer, freeze uncertainty,
   observations and calibration/validation split before fitting, and retain
   independent flat/hill scores. No historical strip result is substituted
   for the new full coupled patch. Quantify numerical error separately from
   physical model error and matching-height/switch-width sensitivity.

All attempts, including failed qualifications, are retained with source,
commands, runtime/dependency fingerprints, allocation and exit status. M8
does not inherit a successful physical gate merely from M7 verification.
