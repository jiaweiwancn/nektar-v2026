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

## Current observations (qualification still in progress)

The first implementation is `3941b83`; the source-bound native campaign uses
`aa854f4c5b05451a81a4e1dd60eaaca98eb547e4` in a separate `m8/source` and
`m8/build`. M7's source/build/runtime evidence and the unrelated Linux main
checkout are preserved. No revised-note or original CLS/AWF source is changed.

- The compiled kernel probe passes 27,907 checks. Independent kernel/source
  qualification passes 37,310 checks at `07363199`; its later stored-report
  replay passes 37,309 (the skipped check is the fresh process exit check).
- The reduced transverse experiment passes 305,254 independent checks. At
  xi=4, p=7/eight radial cells and T=.05 s, the energy ratios are
  .589289214077 for beta=.4, .734610537997 for beta=.45, and
  .333442685761 for beta=.4/GammaScale=.25. This is reference-local reduced
  dynamics; it does not qualify a nonlinear terrain trajectory.
- A fresh directed-interval proof at `11aa94c0` certifies both covariance
  positivity and an incompressible momentum principal lower bound for the
  entire local box with gradient-entry radius .001 and k/epsilon ratios
  [.99,1.01], using GammaScale=1 for both beta choices. GammaScale=3
  jointly certifies radius .01 with those same turbulence bounds. These
  finite boxes are not invariant regions of the coupled PDE.
- The coarse p=3 native projection is outside that box before any step,
  with max normalized gradient error .0642009881 and negative inner
  sampled covariance. The zero-step p=5 error is .00264835 and still has
  negative inner covariance. At p=7 the error is .00010243325, every
  declared finite outer evaluation set fits the GammaScale=1 local box,
  and inner/outer sampled covariance minima are +.08069398/+.17775442.
  These are reference-projection observations, not a positive time map.
- On the broader preregistered box (gradient-entry radius .25, k/epsilon
  ratios [.5,1.5]), exact rational states/directions with outward interval
  stress/tangent evaluation give contradictory *necessary* fixed-GammaScale
  bounds. For beta=.4 they are GammaScale >=2.44640920758 and
  GammaScale <=.177525831862; for beta=.45, >=2.27794509320 and
  <=.156117116287. Thus neither fixed-scalar regularized variant can
  provide both nonnegative principal diffusion and covariance positivity
  at every state of this box. This is stronger than the earlier failure
  of very conservative sufficient-gamma choices, and applies only to the
  declared box/model family. Smaller regimes and changed invariant laws
  remain separate questions.

The approved modified-invariant alternative is therefore being screened
mathematically before any new production operator is adopted. The concrete
screen retains C(q) below qc=3.3, uses the earlier note's C1 power-law tail
with alpha=-qc*C'(qc)/C(qc), and compares constant common beta with a tail
beta(q) that holds beta(q)[q*C(q)]^2 fixed. All coefficient derivatives
must enter the selected stress tangent. Repair of a simple-shear/transverse
factor alone cannot qualify the full tensor.

The first regularized flat manufactured cases pass their actual solver and
independent field/traction checks. Their campaign then exposes a Python
checkpoint-path serialization error; no checkpoint solve was attempted
after that exception. The corrected driver will preserve the original
attempt and replay only successful, hash-bound cases, with every independent
audit repeated; missing checkpoint/restart cases run afresh. Neither this
receipt bug nor an incomplete trajectory is a constitutive diagnosis.
